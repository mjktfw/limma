#  PLOTS

imageplot <- function(z, layout=list(ngrid.r=12,ngrid.c=4,nspot.r=26,nspot.c=26), low=NULL, high=NULL, ncolors=123, zerocenter=NULL, zlim=NULL,...) {
#  Image plot of spotted microarray data
#  Gordon Smyth
#  20 Nov 2001.  Last revised 27 Nov 2001.

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
	old.par <- par(no.readonly = TRUE)
	on.exit(par(old.par))
	par(mar=rep(1,4))
	image(0:(gr*sr),0:(gc*sc),z,zlim=zlim,col=col,axes=FALSE,...)
	for (igrid in 0:gc) lines( c(0,gr*sr), rep(igrid*sc,2) )
	for (igrid in 0:gr) lines( rep(igrid*sr,2), c(0,gc*sc) )
	invisible()
}

plotPrintTipLoess <- function(MA,layout,array=1,span=0.4,...) {
#  MA-plots by print-tip group
#  Gordon Smyth
#  7 April 2003.  Last revised 28 April 2003.

	y <- MA$M[,array]
	x <- MA$A[,array]
	coplot(y~x|factor(gridc(layout))*factor(gridr(layout)),xlab=c("A","Tip Column"),ylab=c("M","Tip Row"),pch=".",span=span,show.given=FALSE,panel=panel.smooth)
	invisible()
}
