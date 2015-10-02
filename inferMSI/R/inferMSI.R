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
                     dimnames = list(names(tumor), names(model$normal)))

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

    # fill in missing data
    scores <- missForest(as.data.frame(scores))$ximp %>%
              as.matrix()


    ######### Infer MSI status #########
    predictions <- cbind(1, scores) %*% model$pred %>%
                   inv.logit()


    ######### Return results #########
    retval <- list(file = file,
                   msi = predictions > cut.off,
                   score = predictions)

    return(retval)
}
