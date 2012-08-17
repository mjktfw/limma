#	LOESS FUNCTIONS

loessFit <- function(y, x, weights=NULL, span=0.3, bin=0.01/(2-is.null(weights)), iterations=4)
#	Fast loess fit for simple x and y
#	This function uses stats:::lowess if no weights and stats:::loess otherwise.
#	It is intended to give a streamlined common interface to the two functions.
#	Gordon Smyth
#	28 June 2003.  Last revised 17 Aug 2012.
{
	n <- length(y)
	out <- list(fitted=rep(NA,n),residuals=rep(NA,n))
	obs <- is.finite(y) & is.finite(x)
	xobs <- x[obs]
	yobs <- y[obs]
	nobs <- length(yobs)
	if(nobs==0) return(out)
	if(is.null(weights)) {
		iter <- iterations-1
		delta = bin * diff(range(xobs))
		lo <- lowess(x=xobs,y=yobs,f=span,iter=iter,delta=delta)
		o <- order(xobs)
		out$fitted[obs][o] <- lo$y
		out$residuals[obs] <- yobs-out$fitted[obs]
	} else {
		weights[!is.finite(weights)] <- 0
		wobs <- weights[obs]
		if(sum(wobs>0) < 4+1/span) {
			fit <- lm.wfit(cbind(1,xobs),yobs,wobs)
		} else {
#			Suppress warning "k-d tree limited by memory"
#			oldopt <- options(warning.expression=expression())
			oldopt <- options(warn=-1)
			on.exit(options(oldopt))
#			fit <- .vsimpleLoess(y=yobs, x=xobs, weights=wobs, span=span, degree=1, cell=bin/span, iterations=iterations)
			fit <- stats:::simpleLoess(y=yobs,x=xobs,weights=wobs,span=span,degree=1,parametric=FALSE,normalize=FALSE,statistics="none",surface="interpolate",cell=bin/span,iterations=iterations,trace.hat="approximate")
		}
		out$fitted[obs] <- fit$fitted
		out$residuals[obs] <- fit$residuals
	}
	out
}

