#	LOESS FUNCTIONS

loessFit <- function(y, x, weights=NULL, span=0.3, bin=NULL, iterations=4, min.weight=1e-5, max.weight=1e5, equal.weights.as.null=TRUE)
#	Fast loess fit for simple x and y
#	This function uses stats:::lowess if no weights and stats:::loess otherwise.
#	It is intended to give a streamlined common interface to the two functions.
#	Gordon Smyth
#	28 June 2003.  Last revised 8 Oct 2012.
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
		if(is.null(bin)) bin <- 0.01
		delta = bin * diff(range(xobs))
		o <- order(xobs)
#		The .C("lowess" call is copied from stats:::lowess
#		lo <- .C("lowess", x = as.double(xobs[o]), as.double(yobs[o]), 
#			nobs, as.double(span), as.integer(iter), as.double(delta), 
#			y = double(nobs), double(nobs), double(nobs), PACKAGE = "stats")
#		out$fitted[obs][o] <- lo$y
#		out$residuals[obs][o] <- lo[[9]]
#		For R 2.16.X
		lo <- lowess(x=xobs,y=yobs,f=span,iter=iter,delta=delta)
		out$fitted[obs][o] <- lo$y
		out$residuals[obs] <- yobs-out$fitted[obs]
	} else {
		if(is.null(bin)) bin <- 0.005
		wobs <- weights[obs]
		wobs[is.na(wobs)] <- 0
		wobs <- pmax(wobs,min.weight)
		wobs <- pmin(wobs,max.weight)
#		Test whether weights are equal
		if(equal.weights.as.null) {
			r <- range(wobs)
			if(r[2]-r[1] < 1e-15) return(Recall(y,x,span=span,bin=bin,iterations=iterations))
		}
#		Count number of observations with positive weights
		if(min.weight>0)
			nwobs <- nobs
		else 
			nwobs <- sum(wobs>0)
		if(nwobs < 4+1/span) {
			if(nwobs>1) {
				fit <- lm.wfit(cbind(1,xobs),yobs,wobs)
			} else {
				fit <- list()
				fit$fitted <- rep(sum(wobs*yobs)/sum(wobs),nobs)
				fit$residuals <- yobs-fit$fitted
			}
		} else {
#			Suppress warning "k-d tree limited by memory"
#			oldopt <- options(warning.expression=expression())
			oldopt <- options(warn=-1)
			on.exit(options(oldopt))
#			fit <- .vsimpleLoess(y=yobs, x=xobs, weights=wobs, span=span, degree=1, cell=bin/span, iterations=iterations)
#			For R 2.16.X
			fit <- stats:::simpleLoess(y=yobs,x=xobs,weights=wobs,span=span,degree=1,parametric=FALSE,normalize=FALSE,statistics="none",surface="interpolate",cell=bin/span,iterations=iterations,trace.hat="approximate")
		}
		out$fitted[obs] <- fit$fitted
		out$residuals[obs] <- fit$residuals
	}
	out
}

