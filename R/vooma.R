vooma <- function(y,design=NULL,correlation,block=NULL,plot=FALSE,span=NULL)
# Linear modelling of microarray data mean-variance modelling at the observational level.
# Creates an EList object for entry to lmFit() etc in the limma pipeline.
# Gordon Smyth and Charity Law
# Created 31 July 2012.  Last modified 5 Aug 2012.
{
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	narrays <- ncol(y)
	ngenes <- nrow(y)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		design <- matrix(1,narrays,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "GrandMean"
	}

#	Fit linear model
	if(is.null(block)) {
		fit <- lm.fit(design,t(y$E))
		mu <- fit$fitted.values
	} else {
		block <- as.vector(block)
		if(length(block)!=narrays) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,narrays,nblocks)==matrix(ub,narrays,nblocks,byrow=TRUE)
		cormatrix <- Z%*%(correlation*t(Z))
		diag(cormatrix) <- 1
		cholV <- chol(cormatrix)
		z <- backsolve(cholV,t(y$E),transpose=TRUE)
		X <- backsolve(cholV,design,transpose=TRUE)
		fit <- lm.fit(X,z)
		mu <- crossprod(cholV,fit$fitted.values)
	}
	s2 <- colMeans(fit$effects[-(1:fit$rank),,drop=FALSE]^2)

#	Fit lowess trend to sqrt-standard-deviations by ave log intensity
	sx <- rowMeans(y$E)
	sy <- sqrt(sqrt(s2))
	if(is.null(span)) if(ngenes<=10) span <- 1 else span <- 0.3+0.7*(10/ngenes)^0.5
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="Average log2 expression",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lines(l,col="red")
	}

#	Make interpolating rule
	f <- approxfun(l, rule=2)

#	Apply trend to individual observations
	w <- 1/f(t(mu))^4
	dim(w) <- dim(y)

#	Output
	y$weights <- w
	y$design <- design
	y$span <- span
	y
}
