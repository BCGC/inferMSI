# trainMSI.R
# Create a model for MSI inference using normal data and cancer samples with known MSI status
# Randy Johnson
# CCR Collaborative Bioinformatics Resrouce at Frederick National Laboratory
# Leidos Biomedical Research, Inc

trainMSI <- function(germline, somatic, msi, markers = NULL, nValidations = 1000, train.prop = 0.5,
                     validate.prop = 0.2, seed = 2308947, verbose = FALSE)
{
    ######### Checks #########
    # ... should add some checks here to be sure variables make sense ... assuming they do for now


    ######### Read .repeatseq files #########
    normal <- lapply(germline, read.repeatseq, markers)
    tumor <- lapply(somatic, read.repeatseq, markers)


    ######### Collect Mixture Distributions #########
    dstn <- list()

    for(i in 1:length(normal[[1]]))
    {
        # gather individual distributions
        tmp <- sapply(normal, `[`, i)

        dstn[[names(normal[[1]])[i]]] <- numeric()

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

    for(i in 1:length(tumor))
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
    scores$msi <- msi

    # just in case...
    scores <- subset(scores, !is.na(msi))


    ######### Training #########

    set.seed(seed)

    # pick training/testing groups
    if(train.prop < 1)
    {
        training <- c(sample((1:dim(scores)[1])[ scores$msi], size = round(sum( scores$msi)*train.prop)),
                      sample((1:dim(scores)[1])[!scores$msi], size = round(sum(!scores$msi)*train.prop)))
        testing <- (1:dim(scores)[1])[-training]
    }else{
        training <- 1:dim(scores)[1]
        testing <- numeric()
    }

    # impute missing scores
    tmp <- missForest(scores[training,names(scores) %in% c(names(dstn), 'msi')])$ximp

    # get response variable column number
    resp <- which(names(tmp) == 'msi')

    tmp[,resp] <- as.numeric(tmp[,resp])
    tmp <- as.matrix(tmp)

    # store cross validation results here
    validate <- matrix(nrow = 1000, ncol = ncol(tmp) + 1)
    colnames(validate) <- c(colnames(tmp), 'score')
    validationScore <- ncol(tmp) + 1

    # run a bunch of models
    for(i in 1:nValidations)
    {
        # sample 80% / 20% for training / testing
        nmsi <- sum(scores$msi[training])
        test <- c(sample((1:dim(tmp)[1])[scores$msi[training]], nmsi*validate.prop),
                  sample((1:dim(tmp)[1])[!scores$msi[training]], (dim(tmp)[1] - nmsi) * validate.prop))
        train <- (1:dim(tmp)[1])[-test]

        # model
        validate[i,-validationScore] <- logitreg(as.numeric(tmp[train, resp]), as.matrix(tmp[train, -resp]))$par

        # validate
        predict <- inv.logit(cbind(1, tmp[test, -resp]) %*% validate[i,-validationScore])

        # calculate likelihood score given correct solution (i.e. the larger this is the worse the prediction)
        validate[i, validationScore] <- sum(ifelse(tmp[test,resp], log10(predict), log10(1 - predict)))
    }


    # pick final model
    model <- list()
    model$pred <- apply(validate[is.finite(validate[,'score']),-validationScore], 2, median)
    model$markers <- markers

    # if verbose show full spectrum of models
    if(verbose)
        if(require(MASS))
            parcoord(validate[,-validationScore], lty = ifelse(!is.finite(validate[,validationScore]), 3, 1),
                     col = ifelse(!is.finite(validate[,validationScore]), 'red', 'black'))


    ######### Testing #########
    if(length(testing) > 0)
    {
        tmp <- as.matrix(missForest(scores[testing,names(scores) %in% c(names(dstn), 'msi')])$ximp)
    }

    # predict (posibly using test data)
    prediction <- inv.logit(cbind(1, tmp[,-resp]) %*% model$pred)

    # collect model metrics via ROC
    roc <- matrix(nrow = length(unique(prediction)) + 2, ncol = 3,
                  dimnames = list(NULL, c('False Positive', 'True Positive', 'Cutoff')))

    roc[,'Cutoff'] <- c(1, unique(prediction[order(-prediction)]), 0)

    for(i in 1:dim(roc)[1])
    {
        roc[i,'False Positive'] <- sum(prediction[tmp[,'msi'] == 0] > roc[i,'Cutoff']) / sum(tmp[,'msi'] == 0)
        roc[i,'True Positive'] <- sum(prediction[tmp[,'msi'] == 1] > roc[i,'Cutoff']) / sum(tmp[,'msi'] == 1)
    }

    model$roc <- roc
    model$auc <- sum((roc[-dim(roc)[1],'Cutoff'] -  roc[-1,'Cutoff']) * roc[-1,'True Positive'])

    if(verbose)
    {
        plot(roc[,1], roc[,2], type = 'l', xlab = 'False Positive', ylab = 'True Positive', main = 'ROC')
        abline(0:1, lty = 2)

        plot(prediction, col = ifelse(tmp[,'msi'] == 1, 'red', 'black'), pch = 20, ylab = 'Prediction')
    }

    return(model)
}

# read repeatseq data from file <f>
read.repeatseq <- function(f, markers = NULL)
{
    if(!file.exists(f))
        stop(paste("File", f, "seems not to exist"))

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
                          retval <- suppressWarnings(as.numeric(sapply(x, `[`, 2)))
                          names(retval) <- sapply(x, `[`, 1)
                          # order by allele size
                          retval <- retval[order(names(retval))]
                          # convert to proportion
                          retval <- retval / sum(retval)
                      })

    ### concordance
    tmp <- strsplit(sapply(tmp, `[`, 2), ' D:')
    concordance <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### total reads
    tmp <- strsplit(sapply(tmp, `[`, 2), ' R:')
    totalReads <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### total reads after filtering
    tmp <- strsplit(sapply(tmp, `[`, 2), ' S:')
    nReads <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### total reads with no CIGAR sequence present
    tmp <- strsplit(sapply(tmp, `[`, 2), ' M:')
    noCIGAR <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### average mapping quality
    tmp <- strsplit(sapply(tmp, `[`, 2), ' GT:')
    qc <- suppressWarnings(as.numeric(sapply(tmp, `[`, 1)))

    ### genotype
    tmp <- strsplit(sapply(tmp, `[`, 2), ' L:')
    geno <- sapply(tmp, `[`, 1)

    ### likelihodd of genotype
    lik <- suppressWarnings(as.numeric(sapply(tmp, `[`, 2)))


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

    if(!is.null(markers))
    {
        for(j in 1:length(retval))
        {
            if(retval[[j]]$siteInfo %in% names(markers))
            {
                names(retval)[j] <- markers[retval[[j]]$siteInfo]
            }else{
                names(retval)[j] <- retval[[j]]$siteInfo
            }
        }
    }else{
        names(retval) <- unlist(sapply(retval, `[`, 'siteInfo'))
    }

    return(retval)
}
