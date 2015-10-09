# inferMSI.R
# Call MSI status given a model from trainMSI
# Randy Johnson
# CCR Collaborative Bioinformatics Resource at Frederick National Laboratory
# Leidos Biomedical Research, Inc


inferMSI <- function(file, model, cut.off)
{
    ######### Checks #########
    # ... should add some checks here to be sure variables make sense ... assuming they do for now


    ######### Read .repeatseq files #########
    repeatseq <- lapply(file, read.repeatseq, model$markers)


    ######### Score samples #########
    scores <- matrix(NA, nrow = length(repeatseq), ncol = length(model$normal),
                     dimnames = list(names(repeatseq), names(model$normal)))

    for(i in 1:length(repeatseq))
    {
        for(j in names(model$normal))
        {
            tmp <- repeatseq[[i]][[j]]$alleles

            matches <- names(tmp) %in% names(model$normal[[j]])
            tmp[matches] <- abs(tmp[matches] - model$normal[[j]][names(tmp)[matches]])
            scores[i,j] <- sum(tmp)
        }
    }


    ######### fill in missing data #########

    # if there are only a few (or no) observed values, fill in NAs with average score among non-msi samples
    nunique <- apply(scores, 2, unique) %>%
               lapply(na.omit) %>%
               sapply(length)

    for(j in which(nunique < 5))
    {
        if(names(nunique)[j] %in% names(model$pred))
            scores[is.na(scores[,j]),j] <- model$meanScore[names(nunique)[j]]
    }

    # among remaining NAs, impute
    scores <- as.matrix(
                        missForest(
                                   as.data.frame(
                                                 scores[,names(model$pred)[-1]]
                                                )
                                  )$ximp
                        )


    ######### Infer MSI status #########
    predictions <- cbind(1, scores[,]) %*% model$pred %>%
                   inv.logit()


    ######### Return results #########
    retval <- data.frame(file = file,
#                         msi = predictions > cut.off,
                         score = predictions)

    return(retval)
}
