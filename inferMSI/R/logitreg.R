# We need this because the Fisher scoring algorithm used by default in
# glm has problems when there is a perfect fit of the data...this happens
# frquently in this application, particularly when data sets are smallish!
# see MASS book pages 197-8 and 445 (4th Edition)
logitreg <- function(y, x, wt = rep(1, length(y)), intercept = TRUE, start = rep(0, p), ...)
{
    # if x has no columns, reutrn odds based only on y
    if(ncol(x) == 0)
        return(mean(y) / (1 - mean(y)))

    fmin <- function(beta, X, y, w)
    {
        p <- plogis(X %*% beta)
        -sum(2 * w * ifelse(y, log(p), log(1-p)))
    }

    gmin <- function(beta, X, y, w)
    {
        eta <- X %*% beta
        p <- plogis(eta)

        -2 * t(w * dlogis(eta) * ifelse(y, 1/p, -1 / (1 - p))) %*% X
    }

    if(is.null(dim(x)))
       dim(x) <- c(length(x), 1)

    dn <- dimnames(x)[[2]]

    if(!length(dn))
        dn <- paste("PC", 1:ncol(x), sep = "")

    p <- ncol(x) + intercept

    if(intercept)
    {
        x <- cbind(1, x)
        dn <- c("(Intercept)", dn)
    }

    if(is.factor(y))
        y <- unclass(y) != 1

    fit <- optim(start, fmin, gmin, X=x, y=y, w=wt, method = "BFGS", hessian = FALSE, ...)

    names(fit$par) <- dn

    invisible(fit)
} # end of logitreg
