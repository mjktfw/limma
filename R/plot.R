##  PLOT.R
##  plot() produces MA plots on expression objects

plot.RGList <- function(x, y, array=1, xlab="A", ylab="M", main=colnames(x)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014
{
	MA <- MA.RG(x[,array])
	plot.MAList(MA,array=1,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend,zero.weights=zero.weights,...)
}

plot.MAList <- function(x, y, array=1, xlab="A", ylab="M", main=colnames(x)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014
{
	MA <- x
	x <- as.matrix(MA$A)[,array]
	y <- as.matrix(MA$M)[,array]
	if(is.null(MA$weights)) w <- NULL else w <- as.matrix(MA$weights)[,array]
	if(missing(status)) status <- MA$genes$Status
	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plot.MArrayLM <- function(x, y, coef=ncol(x), xlab="AveExpr", ylab="logFC", main=colnames(x)[coef], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014
{
	fit <- x
	if(is.null(fit$Amean)) stop("MA-plot not possible because Amean component is absent.")
	x <- fit$Amean
	y <- as.matrix(fit$coef)[,coef]
	if(is.null(fit$weights)) w <- NULL else w <- as.matrix(fit$weights)[,coef]
	if(missing(status)) status <- fit$genes$Status
	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plot.EList <- function(x, y, array=1, xlab="A", ylab="M", main=colnames(x)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth
#	Created 21 March 2014
{
	E <- x
	E$E <- as.matrix(E$E)
	narrays <- ncol(E$E)
	if(narrays < 2) stop("Need at least two arrays")
	if(narrays > 5)
		x <- apply(E$E,1,median,na.rm=TRUE)
	else
		x <- rowMeans(E$E,na.rm=TRUE)
	y <- E$E[,array]-x
	if(is.null(E$weights)) w <- NULL else w <- as.matrix(E$weights)[,array]
	if(missing(status)) status <- E$genes$Status

	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}
