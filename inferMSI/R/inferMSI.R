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
    repeatseq <- lapply(file, read.repeatseq, markers)


    ######### Score samples #########
    scores <- matrix(NA, nrow = length(repeatseq), ncol = length(dstn),
                     dimnames = list(names(tumor), names(dstn)))

    for(i in 1:length(repeatseq))
    {
        for(j in names(dstn))
        {
            tmp <- repeatseq[[i]][[j]]$alleles

            matches <- names(tmp) %in% names(dstn[[j]])
            tmp[matches] <- abs(tmp[matches] - dstn[[j]][names(tmp)[matches]])
            scores[i,j] <- sum(tmp)
        }
    }

    scores <- as.data.frame(scores)

    # fill in missing data
    scores <- missForest(scores[,names(scores) %in% names(dstn)])$ximp
    scores[,!names(scores) %in% names(dstn)] <- 0

    ######### Infer MSI status #########
    if(dim(tmp)[1] == 1)
    {
        predictions <- inv.logit(t(t(c(1, tmp))) %*% model$pred)
    }else{
        predictions <- inv.logit(cbind(1, tmp) %*% model$pred)
    }


    ######### Return results #########
    retval <- list(msi = prediction > cut.off,
                   score = prediction)

    return(retval)
}
