# trainMSI.R
# Create a model for MSI inference using normal data and cancer samples with known MSI status
# Randy Johnson
# CCR Collaborative Bioinformatics Resrouce at Frederick National Laboratory
# Leidos Biomedical Research, Inc

trainMSI <- function(germline, somatic, msi, markers = NULL, verbose = FALSE)
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


    ######### Calculate Microsat Distributions #########

    # drop any markers that didn't sequence properly
    for(i in names(dstn)[which(sapply(dstn, length) == 0)])
        dstn[[i]] <- NULL

    # cumulative density function
    for(i in 1:length(dstn))
    {
        r <- range(as.numeric(names(dstn[[i]])))

        tmp <- rep(as.numeric(names(dstn[[i]])), ceiling(dstn[[i]]*10)) %>%
               densityMclust() %>%
               cdfMclust(seq(r[1] - 0.5, r[2] + 0.5, by = 1))

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

    # drop any that don't have enough unique values
    nunique <- apply(scores[,names(dstn)], 2, unique) %>% # get unique values
               lapply(na.omit) %>%                        # remove NA from unique values
               sapply(length)                             # count number of uniue values

    scores <- scores[,c(names(nunique)[nunique > 4], 'msi')] # want 5 or more unique values

    # impute missing scores
    scores <- missForest(scores[,names(scores) %in% c(names(dstn), 'msi')])$ximp


    ######### Training #########

    # get response variable column number
    resp <- which(names(scores) == 'msi')
    scores <- as.matrix(scores)

    # store cross validation results here
    validate <- matrix(nrow = dim(scores)[1], ncol = ncol(scores) + 2)
    colnames(validate) <- c('(Intercept)', colnames(scores)[-resp], 'score', 'msi')
    validationScore <- (ncol(validate) - 1):ncol(validate) # these columns contain the validationScores

    # train model with leave-one-out validation
    for(i in 1:dim(scores)[1])
    {
        # model
        validate[i,-validationScore] <- logitreg(as.numeric(scores[-i, resp]), as.matrix(scores[-i, -resp]))$par

        # validate
        validate[i,'score'] <- c(1, scores[i, -resp]) %*% validate[i,-validationScore] %>%
                               inv.logit()
        validate[i,'msi'] <- scores[i,resp]
    }

    Pcorrect <- ifelse(validate[,'msi'], validate[,'score'], 1 - validate[,'score'])


    # pick final model
    model <- list()
    model$pred <- validate[is.finite(validate[,'score']),-validationScore] %>% # only keep rows w/ finite scores
                  apply(2, median) # use median to avoid overly influential rows
    model$markers <- markers
    model$normal <- dstn
    model$meanScore <- apply(scores[!msi,names(model$pred)[-1]], 2, mean, na.rm = TRUE)

    # if verbose show full spectrum of models
    if(verbose)
        if(require(MASS))
            parcoord(validate[,-validationScore], col = rgb(1 - Pcorrect, Pcorrect, 0))


    ######### Testing #########

    # predict using final model
    prediction <- inv.logit(cbind(1, scores[,names(model$pred)[-1]]) %*% model$pred)

    # collect model metrics via ROC
    roc <- matrix(nrow = 300, ncol = 3, dimnames = list(NULL, c('False Positive', 'True Positive', 'Cutoff')))

    roc[,'Cutoff'] <- seq(from = 1, to = 0, length = 300)

    for(i in 1:dim(roc)[1])
    {
        roc[i,'False Positive'] <- sum(validate[validate[,'msi'] == 0, 'score'] > roc[i,'Cutoff'])
        roc[i,'True Positive'] <- sum(validate[validate[,'msi'] == 1, 'score'] > roc[i,'Cutoff'])
    }

    roc[,'False Positive'] <- roc[,'False Positive'] / sum(validate[,'msi'] == 0)
    roc[,'True Positive'] <- roc[,'True Positive'] / sum(validate[,'msi'] == 1)

    # only keep those where there is a change
    roc <- rbind(c(0, 0, 1),
                 roc[c(FALSE, roc[-1,'False Positive'] != roc[-300,'False Positive'] |
                              roc[-1,'True Positive'] != roc[-300,'True Positive']),],
                 c(1, 1, 0))

    model$roc <- roc

    if(verbose)
    {
        plot(roc[,1], roc[,2], type = 'l', xlab = 'False Positive', ylab = 'True Positive', main = 'ROC')
        abline(0:1, lty = 2)

        plot(prediction, col = ifelse(scores[,'msi'] == 1, 'red', 'black'), pch = 20, ylab = 'Prediction')
    }

    return(model)
}


scoreMSI <- function(f, somatic, markers = NULL)
{
    ######### Checks #########
    # ... should add some checks here to be sure variables make sense ... assuming they do for now


    ######### Read .repeatseq files #########
    rptsq <- lapply(f, read.repeatseq, markers)
    names(rptsq) <- f


    ######### Collect Mixture Distributions #########
    dstn <- list()

    for(i in 1:length(rptsq[[1]]))
    {
        # gather individual distributions
        tmp <- sapply(rptsq[!somatic], `[`, i)

        dstn[[names(rptsq[[1]])[i]]] <- numeric()

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

    scores <- matrix(NA, nrow = length(rptsq), ncol = length(dstn),
                 dimnames = list(names(rptsq), names(dstn)))

    for(i in 1:length(rptsq))
    {
        for(j in names(dstn))
        {
            tmp <- rptsq[[i]][[j]]$alleles

            matches <- names(tmp) %in% names(dstn[[j]])
            tmp[matches] <- abs(tmp[matches] - dstn[[j]][names(tmp)[matches]])
            scores[i,j] <- sum(tmp)
        }
    }

    scores <- as.data.frame(scores)
}
