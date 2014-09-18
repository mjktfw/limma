##  PLOT.R
##  plot() produces MA plots (aka mean difference plots) on expression objects

plot.RGList <- function(x, y, array=1, xlab="A", ylab="M", main=colnames(x)[array], status=x$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014. Last modified 18 Sep 2014.
{
	MA <- MA.RG(x[,array])
	plot.MAList(x=MA,array=1,xlab=xlab,ylab=ylab,main=main,status=status,zero.weights=zero.weights,...)
}

plot.MAList <- function(x, y, array=1, xlab="A", ylab="M", main=colnames(x)[array], status=x$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014. Last modified 18 Sep 2014.
{
	A <- as.matrix(x$A)[,array]
	M <- as.matrix(x$M)[,array]
	if(!zero.weights && !is.null(x$weights)) {
		w <- as.matrix(x$weights)[,array]
		M[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=A,y=M,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plot.MArrayLM <- function(x, y, coef=ncol(x), xlab="Average log-expression", ylab="log-fold-change", main=colnames(x)[coef], status=x$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014.  Last modified 18 Sep 2014.
{
	if(is.null(x$Amean)) stop("Amean component is absent.")
	logFC <- as.matrix(x$coef)[,coef]
	if(!zero.weights && !is.null(x$weights)) {
		w <- as.matrix(x$weights)[,coef]
		logFC[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=x$Amean,y=logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plot.EList <- function(x, y, array=1, xlab="Average log-expression", ylab="Expression log-ratio (this sample vs others)", main=colnames(x)[array], status=x$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014. Last modified 18 Sep 2014.
{
	E <- as.matrix(E$E)
	if(ncol(E) < 2) stop("Need at least two arrays")

#	Convert array to integer if not already
	j <- 1L:ncol(E)
	names(j) <- colnames(E)
	array <- j[array[1]]

	AveOfOthers <- rowMeans(E[,-array,drop=FALSE],na.rm=TRUE)
	Diff <- E[,array]-AveOfOthers
	Mean <- (E[,array]+AveOfOthers)/2

	if(!zero.weights && !is.null(x$weights)) {
		w <- as.matrix(x$weights)[,array]
		Diff[ is.na(w) | (w <= 0) ] <- NA
	}

	plotWithHighlights(x=Mean,y=Diff,xlab=xlab,ylab=ylab,main=main,status=status,...)
}
