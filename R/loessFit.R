#	LOESS FUNCTIONS

loessFit <- function(y, x, weights=NULL, span=0.3, iterations=4, min.weight=1e-5, max.weight=1e5, equal.weights.as.null=TRUE)
#	Fast lowess fit for univariate x and y allowing for weights
#	Uses lowess() if weights=NULL and locfit.raw() otherwise
#	Gordon Smyth
#	28 June 2003.  Last revised 2 Oct 2013.
{
	n <- length(y)
	out <- list(fitted=rep(NA,n),residuals=rep(NA,n))
	obs <- is.finite(y) & is.finite(x)
	xobs <- x[obs]
	yobs <- y[obs]
	nobs <- length(yobs)
	if(nobs==0) return(out)
	if(is.null(weights)) {
		o <- order(xobs)
		lo <- lowess(x=xobs,y=yobs,f=span,iter=iterations-1)
		out$fitted[obs][o] <- lo$y
		out$residuals[obs] <- yobs-out$fitted[obs]
	} else {
		wobs <- weights[obs]
		wobs[is.na(wobs)] <- 0
		wobs <- pmax(wobs,min.weight)
		wobs <- pmin(wobs,max.weight)
#		Test whether weights are equal
		if(equal.weights.as.null) {
			r <- range(wobs)
			if(r[2]-r[1] < 1e-15) return(Recall(y=y,x=x,span=span,iterations=iterations))
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
#			Check for locfit package
			loaded <- ( "package:locfit" %in% search() )
			if(!loaded) {
				loadresult <- tryCatch(suppressPackageStartupMessages(library("locfit",character.only=TRUE,quietly=TRUE)),error=function(e) e)
				if(inherits(loadresult,"error")) stop("limma:::loessFit with weights requires locfit package, which is not available",call.=FALSE)
			}

#			Weighted lowess with robustifying iterations
		    biweights <- rep(1,nobs)
 			for (i in 1:iterations) {
        		fit <- locfit.raw(x=xobs,y=yobs,weights=wobs*biweights,alpha=span,deg=1)
		        res <- residuals(fit,type="raw")
		        s <- median(abs(res))
		        biweights <- pmax(1-(res/(6*s))^2,0)^2
		    }
		}
		out$fitted[obs] <- fitted(fit)
		out$residuals[obs] <- res
	}
	out
}

