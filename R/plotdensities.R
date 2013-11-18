#  PLOT DENSITIES

plotDensities <- function(object,...)
UseMethod("plotDensities")

plotDensities.RGList <- function(object,log=TRUE,group=NULL,col=NULL,main="RG Densities",...)
#	Plot empirical single-channel densities
#	Original version by Natalie Thorne, 9 September 2003
#	Modified by Gordon Smyth.  Last modified 18 Nov 2013.
{
	object <- backgroundCorrect(object,method="subtract")
	narray <- ncol(object)
	E <- cbind(object$R,object$G)

	if(log) E <- log2(E+1)

	col2 <- col
	if(is.null(group)) {
		group2 <- factor(rep(1:2,c(narray,narray)),labels=c("R","G"))
		if(is.null(col2)) col2 <- c("red","green")
	} else {
		group <- rep(group,narray)
		group2 <- c(group,group)
	}

	plotDensities(E,group=group2,col=col2,main=main)
}

plotDensities.MAList <- function(object,log=TRUE,group=NULL,col=NULL,main="RG Densities",...)
#	Plot empirical single-channel densities
#	Original version by Natalie Thorne, 9 September 2003
#	Modified by Gordon Smyth.  Last modified 18 Nov 2013.
{
	narray <- ncol(object)
	E <- cbind(object$A+object$M/2, object$A-object$M/2)
	if(!log) E <- 2^E

	col2 <- col
	if(is.null(group)) {
		group2 <- factor(rep(1:2,c(narray,narray)),labels=c("R","G"))
		if(is.null(col2)) col2 <- c("red","green")
	} else {
		group <- rep(group,narray)
		group2 <- c(group,group)
	}

	plotDensities(E,group=group2,col=col2,main=main)
}

plotDensities.EListRaw <- function(object,log=TRUE,group=NULL,col=NULL,main=NULL,...)
{
	object <- backgroundCorrect(object,method="subtract")
	E <- object$E
	if(log) E <- log2(E+1)
	plotDensities(E,group=group,col=col,main=main)
}

plotDensities.EList <- function(object,log=TRUE,group=NULL,col=NULL,main=NULL,...)
{
	E <- object$E
	if(!log) E <- 2^E
	plotDensities(E,group=group,col=col,main=main)
}

plotDensities.default <- function(object,group=NULL,col=NULL,main=NULL,...)
#	Plot empirical single-channel densities
#	Gordon Smyth
#	18 Nov 2013.  Last modified 18 Nov 2013.
{
#	Coerce object to matrix
	E <- as.matrix(object)
	narray <- ncol(E)

#	Check group
	if(is.null(group))  group <- colnames(E)
	if(is.null(group))  group <- 1:narray
	group <- as.factor(group)
	ngroup <- nlevels(group)

#	Check col
	if(is.null(col)) col <- 1:ngroup
	col <- rep(col,length=ngroup)

#	Expand cols to number of arrays
	arraycol <- group
	levels(arraycol) <- col
	arraycol <- as.vector(arraycol)

	npoint <- 512
	X <- Y <- matrix(0,npoint,narray)
	for (a in 1:ncol(E)) {
		d <- density(E[,a],n=npoint)
		X[,a] <- d$x
		Y[,a] <- d$y
	}
	matplot(X,Y,xlab="Intensity",ylab="Density",main=main,type="l",col=arraycol,lwd=2,lty=1)
	if(ngroup>1) legend("topleft",lwd=2,legend=levels(group),col=col)
	invisible(list(X=X,Y=Y))
}
