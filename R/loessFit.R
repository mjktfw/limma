#	LOESS FUNCTIONS

loessFit <- function(y, x, weights=NULL, span=0.3, iterations=4L, min.weight=1e-5, max.weight=1e5, equal.weights.as.null=TRUE, method="weightedLowess")
#	Fast lowess fit for univariate x and y allowing for weights
#	Uses lowess() if weights=NULL and weightedLowess() or locfit.raw() otherwise
#	Gordon Smyth
#	28 June 2003.  Last revised 15 Jan 2014.
{
#	Check x and y
	n <- length(y)
	obs <- is.finite(y) & is.finite(x)
	xobs <- x[obs]
	yobs <- y[obs]
	nobs <- length(yobs)

#	If no good obs, exit straight away
	out <- list(fitted=rep(NA,n),residuals=rep(NA,n))
	if(nobs==0) return(out)

#	Check weights
	if(!is.null(weights)) {
		wobs <- weights[obs]
		wobs[is.na(wobs)] <- 0
		wobs <- pmax(wobs,min.weight)
		wobs <- pmin(wobs,max.weight)
#		If weights all equal, treat as NULL
		if(equal.weights.as.null) {
			r <- range(wobs)
			if(r[2]-r[1] < 1e-15) weights <- NULL
		}
	}

	if(is.null(weights)) {
#		No weights, so use classic lowess algorithm
		o <- order(xobs)
		lo <- lowess(x=xobs,y=yobs,f=span,iter=iterations-1)
		out$fitted[obs][o] <- lo$y
		out$residuals[obs] <- yobs-out$fitted[obs]
	} else {
#		Count number of observations with positive weights
		if(min.weight>0)
			nwobs <- nobs
		else 
			nwobs <- sum(wobs>0)
		if(nwobs < 4+1/span) {
#			Two few obs to estimate lowess curve, so use linear regression or weighted average
			if(nwobs>1) {
				fit <- lm.wfit(cbind(1,xobs),yobs,wobs)
				out$fitted[obs] <- fit$fitted
				out$residuals[obs] <- fit$residuals
			} else {
				out$fitted[obs] <- rep(sum(wobs*yobs)/sum(wobs),nobs)
				out$residual[obs] <- yobs-fit$fitted
			}
		} else {
			method <- match.arg(method, c("weightedLowess","locfit","loess"))
			switch(method,
				"weightedLowess" = {
					fit <- weightedLowess(x=xobs,y=yobs,weights=wobs,span=span,iterations=iterations,npts=200)
					out$fitted[obs] <- fit$fitted
					out$residuals[obs] <- fit$residuals
				},
				"locfit" = {
#					Check for locfit package
					loaded <- ( "package:locfit" %in% search() )
					if(!loaded) {
						loadresult <- tryCatch(suppressPackageStartupMessages(library("locfit",character.only=TRUE,quietly=TRUE)),error=function(e) e)
						if(inherits(loadresult,"error")) stop("locfit package not available",call.=FALSE)
					}
#					Weighted lowess with robustifying iterations
				    biweights <- rep(1,nobs)
 					for (i in 1:iterations) {
 		       			fit <- locfit.raw(x=xobs,y=yobs,weights=wobs*biweights,alpha=span,deg=1)
				        res <- residuals(fit,type="raw")
				        s <- median(abs(res))
				        biweights <- pmax(1-(res/(6*s))^2,0)^2
				    }
					out$fitted[obs] <- fitted(fit)
					out$residuals[obs] <- res
				},
				"loess" = {
#					Suppress warning "k-d tree limited by memory"
					oldopt <- options(warn=-1)
					on.exit(options(oldopt))
					bin <- 0.01
					fit <- loess(yobs~xobs,weights=wobs,span=span,degree=1,parametric=FALSE,normalize=FALSE,statistics="approximate",surface="interpolate",cell=bin/span,iterations=iterations,trace.hat="approximate")
					out$fitted[obs] <- fit$fitted
					out$residuals[obs] <- fit$residuals
				}
			)
		}
	}
	out
}

