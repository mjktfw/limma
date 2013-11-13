##  PLOTS-MA.R
##  M-A PLOTS

plotMA <- function(MA,...) UseMethod("plotMA")

plotMA.RGList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
	MA <- MA.RG(MA[,array])
	plotMA(MA,array=1,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend,zero.weights=zero.weights,...)
}

plotMA.MAList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
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

plotMA.MArrayLM <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
	if(is.null(MA$Amean)) stop("MA-plot not possible because Amean component is absent.")
	x <- MA$Amean
	y <- as.matrix(MA$coef)[,array]
	if(is.null(MA$weights)) w <- NULL else w <- as.matrix(MA$weights)[,array]
	if(missing(status)) status <- MA$genes$Status
	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plotMA.EList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
	MA$E <- as.matrix(MA$E)
	narrays <- ncol(MA$E)
	if(narrays < 2) stop("Need at least two arrays")
	if(narrays > 5)
		x <- apply(MA$E,1,median,na.rm=TRUE)
	else
		x <- rowMeans(MA$E,na.rm=TRUE)
	y <- MA$E[,array]-x
	if(is.null(MA$weights)) w <- NULL else w <- as.matrix(MA$weights)[,array]
	if(missing(status)) status <- MA$genes$Status

	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plotMA.default <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
#	Data is assumed to be single-channel
	MA <- as.matrix(MA)
	narrays <- ncol(MA)
	if(narrays < 2) stop("Need at least two arrays")
	if(narrays > 5)
		x <- apply(MA,1,median,na.rm=TRUE)
	else
		x <- rowMeans(MA,na.rm=TRUE)
	y <- MA[,array]-x
	w <- NULL
	if(missing(status)) status <- NULL

	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

.plotMAxy <- function(x, y, xlab="A", ylab="M", main=NULL, xlim=NULL, ylim=NULL, status, values, pch, col, cex, legend=TRUE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
#	Check legend
	legend.position <- "topleft"
	if(!is.logical(legend)) {
		legend.position <- legend
		legend <- TRUE
	}
	legend.position <- match.arg(legend.position,c("bottomright","bottom","bottomleft","left","topleft","top","topright","right","center"))

#	Check xlim and ylim
	if(is.null(xlim)) xlim <- range(x,na.rm=TRUE)
	if(is.null(ylim)) ylim <- range(y,na.rm=TRUE)

#	Setup plot axes
	plot(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,type="n",...)

#	If no status information, just plot points normally
	if(is.null(status) || all(is.na(status))) {
		if(missing(pch)) pch <- 16
		if(missing(cex)) cex <- 0.3
		points(x,y,pch=pch[[1]],cex=cex[1])
		return(invisible())
	}

#	From here, status is not NULL and not all missing

#	Check values
	if(missing(values)) {
		if(is.null(attr(status,"values")))
			values <- names(sort(table(status),decreasing=TRUE))
		else
			values <- attr(status,"values")
	}
	nvalues <- length(values)

#	Plot non-highlighted points
	sel <- !(status %in% values)
	nonhi <- any(sel)
	if(nonhi) points(x[sel],y[sel],pch=16,cex=0.3)

	if(missing(pch)) {
		if(is.null(attr(status,"pch")))
			pch <- rep(16,nvalues)
		else
			pch <- attr(status,"pch")
	}

	if(missing(cex)) {
		if(is.null(attr(status,"cex"))) {
			cex <- rep(1,nvalues)
			if(!nonhi) cex[1] <- 0.3
		} else
			cex <- attr(status,"cex")
	}

	if(missing(col)) {
		if(is.null(attr(status,"col"))) {
			col <- nonhi + 1:nvalues
		} else
			col <- attr(status,"col")
	}

	pch <- rep(pch,length=nvalues)
	col <- rep(col,length=nvalues)
	cex <- rep(cex,length=nvalues)

#	Plot highlighted classes of points
	for (i in 1:nvalues) {
		sel <- status==values[i]
		points(x[sel],y[sel],pch=pch[[i]],cex=cex[i],col=col[i])
	}

	if(legend) {
		if(is.list(pch))
			legend(legend.position,legend=values,fill=col,col=col,cex=0.9)
		else
			legend(legend.position,legend=values,pch=pch,,col=col,cex=0.9)
	}
	invisible()
}

plotMA3by2 <- function(MA, prefix="MA", path=NULL, main=colnames(MA), zero.weights=FALSE, common.lim=TRUE, device="png", ...)
#	Make files of MA-plots, six to a page
#	Gordon Smyth  27 May 2004.  Last modified 9 June 2007.
{
	if(is(MA,"RGList")) MA <- MA.RG(MA)
	if(is.null(path)) path <- "."
	prefix <- file.path(path,prefix)
	narrays <- ncol(MA)
	npages <- ceiling(narrays/6)
	device <- match.arg(device, c("png","jpeg","pdf","postscript"))
	if(device=="png" & !capabilities(what="png")) stop("png not available. Try another device.")
	if(device=="jpeg" & !capabilities(what="jpeg")) stop("jpeg not available. Try another device.")
	fdevice <- get(device)
	if(device=="postscript") ext <- "ps" else ext <- device
	width <- 6.5
	height <- 10
	if(device %in% c("png","jpeg")) {
		width <- width * 140
		height <- height * 140
	}

	if(!zero.weights && !is.null(MA$weights)) MA$M[MA$weights<=0] <- NA
	if(common.lim) {
		xlim <- range(MA$A,na.rm=TRUE)
		ylim <- range(MA$M,na.rm=TRUE)
	} else {
		xlim <- ylim <- NULL
	}
	for (ipage in 1:npages) {
		i1 <- ipage*6-5
		i2 <- min(ipage*6,narrays)
		fdevice(file=paste(prefix,"-",i1,"-",i2,".",ext,sep=""),width=width,height=height)
		par(mfrow=c(3,2))
		for (i in i1:i2) {
			plotMA(MA,array=i,xlim=xlim,ylim=ylim,legend=(i%%6==1),zero.weights=TRUE,main=main[i],...)
		}
		dev.off()
	}
	invisible()
}

plotPrintTipLoess <- function(object,layout,array=1,span=0.4,...) {
#  MA-plots by print-tip group
#  Gordon Smyth
#  7 April 2003.  Last revised 19 March 2004.

	if(is(object,"RGList")) {
		object <- MA.RG(object[,array])
		array <- 1
	}
	if(!is.null(object$printer) && missing(layout)) layout <- object$printer
	y <- object$M[,array]
	x <- object$A[,array]
	if(dev.cur()==1) plot.new()
	df <- data.frame(y=object$M[,array],x=object$A[,array],gr=factor(gridc(layout)),gc=factor(layout$ngrid.r-gridr(layout)+1))
	coplot(y~x|gr*gc,data=na.omit(df),xlab=c("A","Tip Column"),ylab=c("M","Tip Row"),pch=".",span=span,show.given=FALSE,panel=panel.smooth)
}

mdplot <- function(x,...)
#	Mean-different plot
#	Gordon Smyth
#	16 March 2005
{
	d <- x[,1]-x[,2]
	m <- (x[,1]+x[,2])/2
	plot(m,d,xlab="Mean",ylab="Difference",...)
}
