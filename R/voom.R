voom <- function(counts,design=NULL,lib.size=NULL,normalize.method="none",plot=FALSE,...) 
# Linear modelling of count data mean-variance modelling at the observational level.
# Creates an EList object for entry to lmFit() etc in the limma pipeline.
# Gordon Smyth and Charity Law
# Created 22 June 2011.  Last modified 25 Nov 2011.
{
	out <- list()

#	Check inputs
	if(is(counts,"DGEList")) {
		out$genes <- counts$genes
		out$targets <- counts$samples
		if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
		if(is.null(lib.size)) lib.size <- with(counts$samples,lib.size*norm.factors)
		counts <- counts$counts
	} else {
		isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
		if(isExpressionSet) {
			if(length(fData(counts))) out$genes <- fData(counts)
			if(length(pData(counts))) out$targets <- pData(counts)
			counts <- exprs(counts)
		} else {
			counts <- as.matrix(counts)
		}
	}
	if(is.null(design)) {
		design <- matrix(1,ncol(counts),1)
		rownames(design) <- colnames(counts)
		colnames(design) <- "GrandMean"
	}
	if(is.null(lib.size)) lib.size <- colSums(counts)

#	Fit linear model to log2-counts-per-million
	y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
	y <- normalizeBetweenArrays(y,method=normalize.method)
	fit <- lmFit(y,design,...)
	if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)

#	Fit lowess trend to sqrt-standard-deviations by log-count-size
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
	fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
	fitted.logcount <- log2(fitted.count)

#	Apply trend to individual observations
	w <- 1/f(fitted.logcount)^4
	dim(w) <- dim(fitted.logcount)

#	Output
	out$E <- y
	out$weights <- w
	out$design <- design
	out$lib.size <- lib.size
	new("EList",out)
}
