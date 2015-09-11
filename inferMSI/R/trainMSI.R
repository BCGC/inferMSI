# trainMSI.R
# Create a model for MSI inference using normal data and cancer samples with known MSI status
# Randy Johnson
# CCR Collaborative Bioinformatics Resrouce at Frederick National Laboratory
# Leidos Biomedical Research, Inc

trainMSI <- function(germline, somatic, msi, markers = NULL, train.prop = 0.5, validate.prop = 0.2, seed = 2308947)
{
    ######### Checks #########
    # ... should add some checks here to be sure variables make sense ... assuming they do for now


    ######### Read .repeatseq files #########
    normal <- lapply(germline, read.repeatseq, markers)
    tumor <- lapply(somatic, read.repeatseq, markers)


    ######### Collect Mixture Distributions #########
    dstn <- list()

    for(i in 1:length(normal))
    {
        # gather individual distributions
        tmp <- sapply(normal, `[`, i)

        dstn[[normal[[1]][[i]]$siteInfo]] <- numeric()

        # do this by individual
        for(j in 1:length(tmp))
        {
            alleles <- tmp[[j]]$alleles
            alleles <- alleles[!is.na(alleles)]

            if(any(names(alleles) == 'NA'))
                stop()
            if(length(alleles) == 0)
                next

            # weight by the log10 number of reads
            weight <- log10(tmp[[j]]$nReads)

            # calculate distribution
            dstn[[i]][names(alleles)] <- ifelse(is.na(dstn[[i]][names(alleles)]), alleles*weight,
                                                dstn[[i]][names(alleles)] + alleles*weight)
        }
    }


    ######### Calculate CDFs #########

    # drop any markers that didn't sequence properly
    for(i in names(dstn)[which(sapply(dstn, length) == 0)])
        dstn[[i]] <- NULL

    # cumulative density function
    for(i in 1:length(dstn))
    {
        r <- range(as.numeric(names(dstn[[i]])))

        tmp <- cdfMclust(densityMclust(rep(as.numeric(names(dstn[[i]])), ceiling(dstn[[i]]*10))),
                         seq(r[1] - 0.5, r[2] + 0.5, by = 1))

        tmp$x <- with(tmp, x[-length(x)] + 0.5)
        tmp$y <- with(tmp, y[-1] - y[-length(y)])

        dstn[[i]] <- tmp$y
        names(dstn[[i]]) <- tmp$x
    }


    ######### Calculate individual differences among somatic samples #########

    scores <- matrix(NA, nrow = length(tumor), ncol = length(dstn),
                 dimnames = list(names(tumor), names(dstn)))

    for(i in names(tumor))
    {
        for(j in names(dstn))
        {
            tmp <- tumor[[i]][[j]]$alleles

            matches <- names(tmp) %in% names(dstn[[j]])
            tmp[matches] <- abs(tmp[matches] - dstn[[j]][names(tmp)[matches]])
            scores[i,j] <- sum(tmp)
        }
    }

    scores <- as.data.frame(scores)
    scores$pid <- rownames(scores)
    scores <- merge(scores, map, by.x = 'pid', by.y = 'tumor')
    scores <- subset(scores, !is.na(msi))


    ######### Training #########

    set.seed(seed)

    # pick training/testing groups
    if(training.prop < 1)
    {
        training <- c(sample((1:dim(scores)[1])[ msi], size = round(sum( msi)*training.prop))
                      sample((1:dim(scores)[1])[!msi], size = round(sum(!msi)*training.prop)))
        testing <- (1:dim(scores)[1])[-training]
    }else{
        training <- 1:dim(scores)[1]
        testing <- numeric()
    }

    # impute missing scores
    tmp <- missForest(scores[training,names(scores) %in% c(names(dstn), 'msi')])$ximp
    names(tmp)[-length(names(tmp))] <- paste('s', 1:(dim(tmp)[2] - 1), sep = '')

    # get response variable column number
    resp <- which(names(tmp) == 'msi')

    tmp[,resp] <- as.numeric(tmp[,resp])
    tmp <- as.matrix(tmp)

    # store cross validation results here
    validate <- matrix(nrow = 1000, ncol = ncol(tmp) + 1)

    # run a bunch of models
    for(i in 1:1000)
    {
        # sample 80% / 20% for training / testing
        test <- c(sample((1:dim(scores)[1])[scores$msi], 2), sample((1:dim(scores)[1])[!scores$msi], 13))
        train <- (1:dim(scores)[1])[-test]

        # train
        model <- logitreg(as.numeric(tmp[train, resp]), as.matrix(tmp[train, -resp]), hessian = FALSE)

        # test
        predict <- inv.logit(cbind(1, tmp[test, -resp]) %*% model$par)

        # calculate likelihood score given correct solution (i.e. the larger this is the worse the prediction)
        validation.score <- sum(ifelse(tmp[test,resp], log10(predict), log10(1 - predict)))

        # save results
        validate[i,] <- c(model$par, validation.score)
    }


    # pick final model
    model <- apply(validate[is.finite(validate[,11]),], 2, median)


    ######### Testing #########
    if(length(testing) > 0)
    {
    }

    return(model)
}

# read repeatseq data from file <f>
read.repeatseq <- function(f, markers)
{
    # read in main lines of info
    tmp <- readLines(f)
    tmp <- gsub('~', '', tmp[which(substr(tmp, 1, 1) == '~')])

    # break out different parts:
    ### site info
    tmp <- strsplit(tmp, ' REF:')
    siteInfo <- sapply(tmp, `[`, 1)

    ### reference allele
    tmp <- strsplit(sapply(tmp, `[`, 2), ' A:')
    ref <- sapply(tmp, `[`, 1)

    ### observed alleles
    tmp <- strsplit(sapply(tmp, `[`, 2), ' C:')
    obs <- sapply(tmp, `[`, 1)
    alleles <- lapply(strsplit(sapply(tmp, `[`, 1), '] ', fixed = TRUE), function(x)
                      {
                          # split observations apart
                          x <- strsplit(gsub(']', '', x, fixed = TRUE), '[', fixed = TRUE)
                          retval <- as.numeric(sapply(x, `[`, 2))
                          names(retval) <- sapply(x, `[`, 1)
                          # order by allele size
                          retval <- retval[order(names(retval))]
                          # convert to proportion
                          retval <- retval / sum(retval)
                      })

    ### concordance
    tmp <- strsplit(sapply(tmp, `[`, 2), ' D:')
    concordance <- as.numeric(sapply(tmp, `[`, 1))

    ### total reads
    tmp <- strsplit(sapply(tmp, `[`, 2), ' R:')
    totalReads <- as.numeric(sapply(tmp, `[`, 1))

    ### total reads after filtering
    tmp <- strsplit(sapply(tmp, `[`, 2), ' S:')
    nReads <- as.numeric(sapply(tmp, `[`, 1))

    ### total reads with no CIGAR sequence present
    tmp <- strsplit(sapply(tmp, `[`, 2), ' M:')
    noCIGAR <- as.numeric(sapply(tmp, `[`, 1))

    ### average mapping quality
    tmp <- strsplit(sapply(tmp, `[`, 2), ' GT:')
    qc <- as.numeric(sapply(tmp, `[`, 1))

    ### genotype
    tmp <- strsplit(sapply(tmp, `[`, 2), ' L:')
    geno <- sapply(tmp, `[`, 1)

    ### likelihodd of genotype
    lik <- as.numeric(sapply(tmp, `[`, 2))


    ##### put everything into a list #####
    retval <- list()
    for(j in 1:length(siteInfo))
    {
        retval[[j]] <- list(siteInfo = siteInfo[j],
                            ref = ref[j],
                            alleles = alleles[[j]],
                            concordance = concordance[j],
                            totalReads = totalReads[j],
                            nReads = nReads[j],
                            noCIGAR = noCIGAR[j],
                            qc = qc[j],
                            geno = geno[j],
                            lik = lik[j])
    }

    names(retval) <- unlist(sapply(retval, `[`, 'siteInfo'))

    return(retval)
}
