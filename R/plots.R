#  PLOTS

imageplot <- function(z, layout=list(ngrid.r=12,ngrid.c=4,nspot.r=26,nspot.c=26), low=NULL, high=NULL, ncolors=123, zerocenter=NULL, zlim=NULL, mar=rep(1,4), ...) {
#  Image plot of spotted microarray data
#  Gordon Smyth
#  20 Nov 2001.  Last revised 8 Feb 2004.

#  Check input
	gr <- layout$ngrid.r
	gc <- layout$ngrid.c
	sr <- layout$nspot.r
	sc <- layout$nspot.c
	if(is.null(gr)||is.null(gc)||is.null(sr)||is.null(sc)) stop("Layout needs to contain components ngrid.r, ngrid.c, nspot.r and spot.c")
	if(length(z) != gr*gc*sr*sc) stop("Number of image spots does not agree with layout dimensions")

#  Check colours
	if(is.character(low)) low <- col2rgb(low)/255
	if(is.character(high)) high <- col2rgb(high)/255
	if(!is.null(low) && is.null(high)) high <- c(1,1,1) - low
	if(is.null(low) && !is.null(high)) low <- c(1,1,1) - high

#  Is zlim preset?
	if(!is.null(zlim)) {
		z <- pmax(zlim[1],z)
		z <- pmin(zlim[2],z)
	}

#  Plot differential expression from "green" to "red" or plot one variable from "white" to "blue"?
	zr <- range(z,na.rm=TRUE)
	zmax <- max(abs(zr))
	zmin <- zr[1]
	if(is.null(zerocenter)) zerocenter <- (zmin < 0)
	if(zerocenter) {
		if(is.null(low)) low <- c(0,1,0)
		if(is.null(high)) high <- c(1,0,0)
		if(is.null(zlim)) zlim <- c(-zmax,zmax)
	} else {
		if(is.null(low)) low <- c(1,1,1)
		if(is.null(high)) high <- c(0,0,1)
		if(is.null(zlim)) zlim <- c(zmin,zmax)
	}

#  Now make the plot
	col <- rgb( seq(low[1],high[1],len=ncolors), seq(low[2],high[2],len=ncolors), seq(low[3],high[3],len=ncolors) )
	dim(z) <- c(sc,sr,gc,gr)
	z <- aperm(z,perm=c(2,4,1,3))
	dim(z) <- c(gr*sr,gc*sc)
	old.par <- par(mar=mar)
	on.exit(par(old.par))
	image(0:(gr*sr),0:(gc*sc),z,zlim=zlim,col=col,axes=FALSE,...)
	for (igrid in 0:gc) lines( c(0,gr*sr), rep(igrid*sc,2) )
	for (igrid in 0:gr) lines( rep(igrid*sr,2), c(0,gc*sc) )
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
