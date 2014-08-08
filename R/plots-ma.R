##  PLOTS-MA.R
##  M-A PLOTS

plotMA <- function(MA,...) UseMethod("plotMA")

plotMA.RGList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status=NULL, values=NULL, pch=NULL, col=NULL, cex=NULL, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
	MA <- MA.RG(MA[,array])
	plotMA(MA,array=1,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend,zero.weights=zero.weights,...)
}

plotMA.MAList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status=NULL, values=NULL, pch=NULL, col=NULL, cex=NULL, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 23 April 2013.
{
	x <- as.matrix(MA$A)[,array]
	y <- as.matrix(MA$M)[,array]
	if(is.null(MA$weights)) w <- NULL else w <- as.matrix(MA$weights)[,array]
	if(is.null(status)) status <- MA$genes$Status
	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plotMA.MArrayLM <- function(MA, coef=ncol(MA), xlab="AveExpr", ylab="logFC", main=colnames(MA)[coef], xlim=NULL, ylim=NULL, status=NULL, values=NULL, pch=NULL, col=NULL, cex=NULL, legend=TRUE, zero.weights=FALSE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 21 March 2014.
{
	if(is.null(MA$Amean)) stop("MA-plot not possible because Amean component is absent.")
	x <- MA$Amean
	y <- as.matrix(MA$coef)[,coef]
	if(is.null(MA$weights)) w <- NULL else w <- as.matrix(MA$weights)[,coef]
	if(is.null(status)) status <- MA$genes$Status
	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plotMA.EList <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status=NULL, values=NULL, pch=NULL, col=NULL, cex=NULL, legend=TRUE, zero.weights=FALSE, ...)
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
	if(is.null(status)) status <- MA$genes$Status

	if(!is.null(w) && !zero.weights) {
		i <- is.na(w) | (w <= 0)
		y[i] <- NA
	}
	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

plotMA.default <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA)[array], xlim=NULL, ylim=NULL, status=NULL, values=NULL, pch=NULL, col=NULL, cex=NULL, legend=TRUE, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 8 August 2014.
{
#	Data is assumed to be single-channel
	MA <- as.matrix(MA)
	narrays <- ncol(MA)
	if(narrays<2) stop("Need at least two columns")
	array <- as.integer(array[1L])
	Ave <- rowMeans(MA[,-array,drop=FALSE],na.rm=TRUE)
	y <- MA[,array]-Ave
	x <- (MA[,array]+Ave)/2

	.plotMAxy(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,status=status,values=values,pch=pch,col=col,cex=cex,legend=legend, ...)
}

# Call this plotWithHighlights and document?

.plotMAxy <- function(x, y, xlab="A", ylab="M", main=NULL, xlim=NULL, ylim=NULL, status=NULL, values=NULL, pch=NULL, col=NULL, cex=NULL, legend=TRUE, pch0=16, col0="black", cex0=0.3, ...)
#	MA-plot with color coding for controls
#	Gordon Smyth 7 April 2003, James Wettenhall 27 June 2003.
#	Last modified 13 April 2014.
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
		points(x,y,pch=pch0,cex=cex0)
		return(invisible())
	}

#	From here, status is not NULL and not all NA

#	Check values
#	Default is to set the most frequent status value as background, and to highlight all other status values in order of frequency
	if(is.null(values)) values <- attr(status,"values")
	if(is.null(values)) {
		status.values <- names(sort(table(status),decreasing=TRUE))
		status <- as.character(status)
		values <- status.values[-1]
	}
	nvalues <- length(values)
	if(nvalues==0L) {
		points(x,y,pch=pch0,cex=cex0)
		return(invisible())
	}

#	From here, values has positive length

#	Plot non-highlighted points
	bg <- !(status %in% values)
	nonhi <- any(bg)
	if(nonhi) points(x[bg],y[bg],pch=pch0,cex=cex0)

#	Check parameters for plotting highlighted points

	if(is.null(pch)) pch <- attr(status,"pch")
	if(is.null(pch)) pch <- pch0
	pch <- rep(pch,length=nvalues)

	if(is.null(cex)) cex <- attr(status,"cex")
	if(is.null(cex)) cex <- 1
	cex <- rep(cex,length=nvalues)

	if(is.null(col)) col <- attr(status,"col")
	if(is.null(col)) col <- nonhi + 1L:nvalues
	col <- rep(col,length=nvalues)

#	Plot highlighted points
	for (i in 1:nvalues) {
		sel <- status==values[i]
		points(x[sel],y[sel],pch=pch[[i]],cex=cex[i],col=col[i])
	}

	if(legend) {
		if(nonhi) {
#			Include background value in legend
			bg.value <- unique(status[bg])
			if(length(bg.value) > 1) bg.value <- "Other"
			values <- c(bg.value,values)
			pch <- c(pch0,pch)
			col <- c(col0,col)
			cex <- c(cex0,cex)
		}
		h <- cex>0.5
		cex[h] <- 0.5+0.8*(cex[h]-0.5)
		if(is.list(pch))
			legend(legend.position,legend=values,fill=col,col=col,cex=0.9,pt.cex=cex)
		else
			legend(legend.position,legend=values,pch=pch,,col=col,cex=0.9,pt.cex=cex)
	}
	invisible()
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
