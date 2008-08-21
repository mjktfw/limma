###  treat.R

treat <- function(fit, lfc=0)
#  Moderated t-statistics with threshold
#  Davis McCarthy, Gordon Smyth
#  25 July 2008
{
    coefficients <- as.matrix(fit$coefficients)
    stdev.unscaled <- as.matrix(fit$stdev.unscaled)
    sigma <- fit$sigma
    df.residual <- fit$df.residual
    if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
        is.null(df.residual)) 
        stop("No data, or argument is not a valid lmFit object")
    if (all(df.residual == 0)) 
        stop("No residual degrees of freedom in linear model fits")
    if (all(!is.finite(sigma))) 
        stop("No finite residual standard deviations")
    sv <- squeezeVar(sigma^2, df.residual)
    fit$s2.prior <- sv$var.prior
    fit$s2.post <- sv$var.post
    df.total <- df.residual + sv$df.prior
    fc.up <- (coefficients > lfc)
    fc.down <- (coefficients < -lfc)
    acoef <- abs(coefficients)
    se <- stdev.unscaled*sqrt(fit$s2.post)
    tstat.right <- (acoef-lfc)/se
    tstat.left <- (acoef+lfc)/se
#    fit$t <- array(0,dim(coefficients),dimnames=dimnames(coefficients))
#    fit$t[fc.up] <- tstat.right[fc.up]
#    fit$t[fc.down] <- -tstat.right[fc.down]
    fit$p.value <- pt(tstat.right, df = df.total,lower=FALSE) + pt(tstat.left,df = df.total,lower=FALSE)
    fit
}
