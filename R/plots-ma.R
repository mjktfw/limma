##  PLOTS-MA.R
##  M-A PLOTS

plotMA <- function(MA,...) UseMethod("plotMA")

plotMA.RGList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], status=MA$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 18 Sep 2014.
{
	MA <- MA.RG(MA[,array])
	plotMA.MAList(MA=MA,array=1,xlab=xlab,ylab=ylab,main=main,status=status,zero.weights=zero.weights,...)
}

plotMA.MAList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], status=MA$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 18 Sep 2014.
{
	A <- as.matrix(MA$A)[,array]
	M <- as.matrix(MA$M)[,array]
	if(!zero.weights && !is.null(MA$weights)) {
		w <- as.matrix(MA$weights)[,array]
		M[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=A,y=M,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMA.MArrayLM <- function(MA, coef=ncol(MA), xlab="Average log-expression", ylab="log-fold-change", main=colnames(MA)[coef], status=MA$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 18 Sep 2014.
{
	if(is.null(MA$Amean)) stop("Amean component is absent.")
	logFC <- as.matrix(MA$coef)[,coef]
	if(!zero.weights && !is.null(MA$weights)) {
		w <- as.matrix(MA$weights)[,array]
		logFC[ is.na(w) | (w <= 0) ] <- NA
	}
	plotWithHighlights(x=MA$Amean,y=logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMA.EList <- function(MA, array=1, xlab="Average log-expression", ylab="Expression log-ratio (this sample vs others)", main=colnames(MA)[array], status=MA$genes$Status, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 18 Sep 2014.
{
	E <- as.matrix(MA$E)
	if(ncol(E) < 2) stop("Need at least two arrays")

#	Convert array to integer if not already
	j <- 1L:ncol(E)
	names(j) <- colnames(E)
	array <- j[array[1]]

	AveOfOthers <- rowMeans(E[,-array,drop=FALSE],na.rm=TRUE)
	Diff <- E[,array]-AveOfOthers
	Mean <- (E[,array]+AveOfOthers)/2

	if(!zero.weights && !is.null(MA$weights)) {
		w <- as.matrix(MA$weights)[,array]
		Diff[ is.na(w) | (w <= 0) ] <- NA
	}

	plotWithHighlights(x=Mean,y=Diff,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMA.default <- function(MA, array=1, xlab="Average log-expression", ylab="Expression log-ratio (this sample vs others)", main=colnames(MA)[array], status=NULL, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 18 Sep 2014.
{
#	Data is assumed to be single-channel
	MA <- as.matrix(MA)
	narrays <- ncol(MA)
	if(narrays<2) stop("Need at least two columns")
	array <- as.integer(array[1L])
	Ave <- rowMeans(MA[,-array,drop=FALSE],na.rm=TRUE)
	y <- MA[,array]-Ave
	x <- (MA[,array]+Ave)/2

	plotWithHighlights(x,y,xlab=xlab,ylab=ylab,main=main,status=status, ...)
}

plotMA3by2 <- function(MA, prefix="MA", path=NULL, main=colnames(MA), zero.weights=FALSE, common.lim=TRUE, device="png", ...)
#	Make files of MA-plots, six to a page
#	Gordon Smyth  27 May 2004.  Last modified 12 Feb 2014.
{
	if(is(MA,"RGList")) MA <- MA.RG(MA)
	if(is(MA,"EListRaw")) {
		A <- rowMeans(log2(MA$E))
		MA$A <- array(A,dim(MA))
		MA$M <- log2(MA$E)-MA$A
		MA <- new("MAList",unclass(MA))
	}
	if(is(MA,"EList")) {
		A <- rowMeans(MA$E)
		MA$A <- array(A,dim(MA))
		MA$M <- MA$E-MA$A
		MA <- new("MAList",unclass(MA))
	}
	if(is.matrix(MA)) {
		A <- rowMeans(MA)
		MA <- new("MAList",list(M=MA-A, A=array(A,dim(MA))))
	}
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
#	Mean-difference plot
#	Gordon Smyth
#	16 March 2005
{
	d <- x[,1]-x[,2]
	m <- (x[,1]+x[,2])/2
	plot(m,d,xlab="Mean",ylab="Difference",...)
}
