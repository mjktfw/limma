voom <- function(counts,design=NULL,lib.size=NULL,normalize.method="none",plot=FALSE,...) 
# Linear modelling of count data mean-variance modelling at the observational level.
# Creates an EList object for entry to lmFit() etc in the limma pipeline.
# Gordon Smyth and Charity Law
# Created 22 June 2011.  Last modified 25 Sep 2011.
{
#	Fit linear model to log2-counts-per-million
	counts <- as.matrix(counts)
	if(is.null(design)) {
		design <- matrix(1,ncol(counts),1)
		rownames(design) <- colnames(counts)
		colnames(design) <- "GrandMean"
	}
	if(is.null(lib.size)) lib.size <- colSums(counts)
	y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
	y <- normalizeBetweenArrays(y,method=normalize.method)
	fit <- lmFit(y,design,...)
	if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)

#	Fit lowess trend to sqrt-standard-deviations by log-count-size
#	sx <- 2^(0.25*(fit$Amean+mean(log2(lib.size+0.5))-log2(1e6)))
	sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)
	sy <- sqrt(fit$sigma)
	allzero <- rowSums(counts)==0
	if(any(allzero)) {
		sx <- sx[!allzero]
		sy <- sy[!allzero]
	}
	l <- lowess(sx,sy,f=0.5)
	if(plot) {
		plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lines(l,col="red")
	}

#	Make interpolating rule
#	Special treatment of zero counts is now removed;
#	instead zero counts get same variance as smallest gene average.
#	l$x <- c(0.5^0.25, l$x)
#	l$x <- c(log2(0.5), l$x)
#	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
#	var0 <- max(var0,1e-6)
#	l$y <- c(var0, l$y)
	f <- approxfun(l, rule=2)

#	Find individual quarterroot fitted counts
	fitted.cpm <- 2^(fit$coef %*% t(fit$design))
	fitted.c <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
#	fitted.qrc <- fitted.c^0.25
	fitted.qrc <- log2(fitted.c)

#	Apply trend to individual observations
	w <- 1/f(fitted.qrc)^4
	dim(w) <- dim(fitted.qrc)
	new("EList",list(E=y,weights=w,design=design,lib.size=lib.size))
}

