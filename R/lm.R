#  LINEAR MODELS

lmFit <- function(object,design=NULL,contrasts=NULL,ndups=1,spacing=1,correlation=0.75,weights=NULL,method="ls",...) {
	gene <-
	if(is(object,"MAList")) {
		if(!is.null(object$design)) design <- object$design
		if(!is.null(object$contrasts)) contrasts <- object$contrasts
		if(!is.null(object$printerlayout)) {
			ndups <- printerlayout$ndups
			spacing <- printerlayout$spacing
		}
		weights <- object$weights
		if(!is.null(gal)) gene <- as.character(gal[,5])
		object <- as.matrix(object$M)
	}
	if(is(object,"marrayNorm")) {
#		don't use accessor function so don't have to require marrayClasses
		weights <- object@maW
		if(length(w) == 0) weights <- NULL
		object <- object@maM
	}
	if(is(object,"exprSet")) {
#		don't use accessor function so don't have to require Biobase
  		weights <- object@se.exprs
		if(length(weights) == 0)
			weights <- NULL
		else
			weights <- 1/pmax(weights,1e-14)^2
		object <- object@exprs
	}
	if(is.null(design)) design <- matrix(1,ncol(object),1)
	if(is.na(correlation) || is.null(correlation)) correlation <- numeric(0)
	method <- match.arg(method,c("ls","robust"))
	if(method=="robust")
		fit <- rlm.series(object,design=design,ndups=ndups,spacing=spacing,weights=weights,...)
	else
		if(ndups < 2 || correlation==0 || length(correlation)==0)
			fit <- lm.series(object,design=design,ndups=ndups,spacing=spacing,weights=weights,...)
		else
			fit <- gls.series(object,design=design,ndups=ndups,spacing=spacing,correlation=correlation,weights=weights,...)
	if(is.null(contrasts))
		contrasts <- matrix(0,0,0)
	else
		fit <- contrasts.fit(fit,contrasts)	
	eb <- ebayes(fit)
	new("MArrayLM",
		design=design,
		contrasts=as.matrix(contrasts),
		coefficients=as.matrix(fit$coefficients),
		stdev.unscaled=as.matrix(fit$stdev.unscaled),
		s2.residual=fit$sigma^2,
		df.residual=fit$df.residual,
		s2.prior=eb$s2.prior,
		df.prior=eb$df.prior,
		s2.post=eb$s2.post,
		tstat=as.matrix(eb$t),
		varcoef.prior=eb$var.prior,
		call=match.call()
	)
}

unwrapdups <- function(M,ndups=2,spacing=1) {
#	Unwrap M matrix for a series of experiments so that all spots for a given gene are in one row
#	Gordon Smyth
#	18 Jan 2002. Last revised 2 Nov 2002.

	if(ndups==1) return(M)
	M <- as.matrix(M)
	nspots <- dim(M)[1]
	nslides <- dim(M)[2]
	ngroups <- nspots / ndups / spacing
	dim(M) <- c(spacing,ndups,ngroups,nslides)
	M <- aperm(M,perm=c(1,3,2,4))
	dim(M) <- c(spacing*ngroups,ndups*nslides)
	M
}

uniquegenelist <- function(genelist,ndups=2,spacing=1) {
#	Eliminate entries in genelist for duplicate spots
#	Gordon Smyth
#	2 Nov 2002.  Last revised 12 Apr 2003

	if(ndups <= 1) return(genelist)
	if(is.null(dim(genelist))) dim(genelist) <- c(length(genelist),1)
	index <- drop(unwrapdups(1:nrow(genelist),ndups=ndups,spacing=spacing)[,1])
	drop(genelist[index,])
}

lm.series <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL)
{
#	Fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	18 Apr 2002. Revised 18 Jan 2003.

	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	if(is.vector(design)) dim(design) <- c(narrays,1)
	nbeta <- ncol(design)
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(NULL,colnames(design)))
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		if(sum(obs) > 0) {
			X <- design[obs,,drop=FALSE]
			y <- y[obs]
			if(is.null(weights))
				out <- lm.fit(X,y)
			else {
				w <- as.vector(weights[i,obs])
				out <- lm.wfit(X,y,w)
			}
			est <- !is.na(out$coef)
			beta[i,] <- out$coef
			stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
			df.residual[i] <- out$df.residual
			if(df.residual[i] > 0)
				if(is.null(weights))
					sigma[i] <- sqrt(sum(out$residuals^2)/out$df.residual)
				else
					sigma[i] <- sqrt(sum(w*out$residuals^2)/out$df.residual)
		}
	}
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual)
}

rlm.series <- function(M,design=NULL,ndups=1,spacing=spacing,weights=NULL,...)
{
#	Robustly fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	20 Mar 2002.  Last revised 2 Nov 2002.

	require(MASS) # need rlm.default
	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	if(is.null(dim(design))) dim(design) <- c(narrays,1)
	nbeta <- ncol(design)
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta)
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		X <- design[obs,,drop=FALSE]
		y <- y[obs]
		if(is.null(weights))
			w <- rep(1,length(y))
		else
			w <- as.vector(weights[i,obs])
		if(length(y) > nbeta) {
			out <- rlm.default(x=X,y=y,weights=w,...)
			beta[i,] <- coef(out)
			stdev.unscaled[i,] <- sqrt(diag(chol2inv(out$qr$qr)))
			df.residual[i] <- length(y) - out$rank
			if(df.residual[i] > 0) sigma[i] <- out$s
		}
	}
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual)
}

gls.series <- function(M,design=NULL,ndups=2,spacing=1,correlation=NULL,weights=NULL,...)
{
#	Fit linear model for each gene to a series of microarrays.
#	Fit is by generalized least squares allowing for correlation between duplicate spots.
#	Gordon Smyth
#	11 May 2002.  Last revised 20 Feb 2003.

	if(ndups<2) {
		warning("No duplicates: correlation between duplicates set to zero")
		ndups <- 1
		correlation <- 0
	}
	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	if(is.null(dim(design))) dim(design) <- c(narrays,1)
	if(nrow(design) != narrays) stop("Number of rows of design matrix does not match number of arrays")
	if(is.null(correlation)) correlation <- dupcor.series(M,design,ndups,...)$cor
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		M[weights < 1e-15 ] <- NA
		weights[weights < 1e-15] <- NA
	}
	nbeta <- ncol(design)
	coef.names <- colnames(design)
	M <- unwrapdups(M,ndups=ndups,spacing=spacing)
	ngenes <- nrow(M)
	if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	design <- design %x% rep(1,ndups)
	cormatrix <- diag(rep(correlation,narrays)) %x% array(1,c(ndups,ndups))
	diag(cormatrix) <- 1
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(NULL,coef.names))
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		y <- y[o]
		n <- length(y)
		if(n > 0) {
			X <- design[o,,drop=FALSE]
			V <- cormatrix[o,o]
			if(!is.null(weights)) {
				wrs <- 1/sqrt(drop(weights[i,o]))
				V <- wrs * t(wrs * t(V))
			}
			cholV <- chol(V)
			y <- backsolve(cholV,y,transpose=TRUE)
			if(all(X==0)) {
				df.residual[i] <- n
				sigma[i] <- array(1/n,c(1,n)) %*% y^2
			} else {
				X <- backsolve(cholV,X,transpose=TRUE)
				out <- lm.fit(X,y)
				est <- !is.na(out$coefficients)
				beta[i,] <- out$coefficients
				stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
				df.residual[i] <- out$df.residual
				if(df.residual[i] > 0)
					sigma[i] <- sqrt( array(1/out$df.residual,c(1,n)) %*% out$residuals^2 )
			}
		}
	}
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual,correlation=correlation)
}

dupcor.series <- function(M,design=rep(1,ncol(M)),ndups=2,spacing=1,initial=0.8,trim=0.15,weights=NULL)
{
#	Estimate the correlation between duplicates given a series of arrays
#	Gordon Smyth
#	25 Apr 2002. Last revised 28 Jan 2003.

	M <- as.matrix(M)
	if(ndups<2) {
		warning("No duplicates: correlation between duplicates not estimable")
		return( list(cor=NA,cor.genes=rep(NA,nrow(M))) )
	}
	require( "nlme" ) # need gls function
	narrays <- ncol(M)
	if(is.vector(design)) dim(design) <- c(narrays,1)
	if(nrow(design) != narrays) stop("Number of rows of design matrix does not match number of arrays")
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		M[weights < 1e-15 ] <- NA
		weights[weights < 1e-15] <- NA
	}
	nbeta <- ncol(design)
	M <- unwrapdups(M,ndups=ndups,spacing=spacing)
	ngenes <- nrow(M)
	if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	design <- design %x% rep(1,ndups)
	Array <- rep(1:narrays,rep(ndups,narrays))
	rho <- rep(NA,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		if(any(!is.finite(y))) y[!is.finite(y)] <- NA
		if(any(diff(Array[is.finite(y)])==0) && sum(!is.na(y)) > nbeta+1)
		if(!is.null(weights)) {
			w <- 1/drop(weights[i,])
			rho[i] <- coef(gls(y~design-1,correlation=corCompSymm(initial,form=~1|Array,fixed=FALSE),weights=~w,na.action=na.omit,control=list(singular.ok=TRUE,returnObject=TRUE,apVar=FALSE))$modelStruct,FALSE)
		} else
			rho[i] <- coef(gls(y~design-1,correlation=corCompSymm(initial,form=~1|Array,fixed=FALSE),na.action=na.omit,control=list(singular.ok=TRUE,returnObject=TRUE,apVar=FALSE))$modelStruct,FALSE)
	}
	rhom <- tanh(mean(atanh(rho),trim=trim,na.rm=TRUE))
	list(cor=rhom,cor.genes=rho)
}

contrasts.fit <- function(fit,contrasts) {
#	Extract contrast information from oneway linear model fit
#	Gordon Smyth
#	13 Oct 2002

	out <- fit
	out$coefficients <- fit$coefficients %*% contrasts
	out$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
	out
}
