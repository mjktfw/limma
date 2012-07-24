##  PLOTMDS.R

#	Class to hold multidimensional scaling output
setClass("MDS",representation("list"))

setMethod("show","MDS",function(object) {
	cat("An object of class MDS\n")
	print(unclass(object))
})

plotMDS <- function(x,...) UseMethod("plotMDS")

plotMDS.MDS <- function(x,labels=colnames(x$distance.matrix),col=NULL,cex=1,dim.plot=x$dim.plot,xlab=paste("Dimension",dim.plot[1]),ylab=paste("Dimension",dim.plot[2]),...)
#	Method for MDS objects
#	Create a new plot using MDS coordinates or distances previously created
#	Gordon Smyth
#	21 May 2011.  Last modified 24 June 2012.
{
#	Are new dimensions requested?
	if(!all(dim.plot==x$dim.plot)) {
		ndim <- max(dim.plot)
		if(ndim > ncol(x$cmdscale.out)) x$cmdscale.out <- cmdscale(as.dist(x$distance.matrix),k=ndim)
		x$x <- x$cmdscale.out[,dim.plot[1]]
		x$y <- x$cmdscale.out[,dim.plot[2]]
	}

#	Estimate width of labels in plot coordinates.
#	Estimate will be ok for default plot width, but maybe too small for smaller plots.
	if(is.null(labels)) labels <- 1:length(x$x)
	StringRadius <- 0.01*cex*nchar(labels)
	left.x <- x$x-StringRadius
	right.x <- x$x+StringRadius

#	Redo plot
	plot(c(left.x,right.x),c(x$y,x$y),type="n",xlab=xlab,ylab=ylab,...)
	text(x$x,x$y,labels=labels,col=col,cex=cex)
	return(invisible(x))
}

plotMDS.default <- function(x,top=500,labels=colnames(x),col=NULL,cex=1,dim.plot=c(1,2),ndim=max(dim.plot),gene.selection="pairwise",xlab=paste("Dimension",dim.plot[1]),ylab=paste("Dimension",dim.plot[2]),...)
#	Multi-dimensional scaling with top-distance
#	Di Wu and Gordon Smyth
#	19 March 2009.  Last modified 8 July 2012.
{
	x <- as.matrix(x)

#	Remove rows with missing or Inf values
	ok <- is.finite(x)
	if(!all(ok)) x <- x[apply(ok,1,all),]

	if(is.null(labels)) labels<-1:dim(x)[2]

	nprobes <- nrow(x)
	nsamples <- ncol(x)
	if(ndim < 2) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("Two few samples")
	top <- min(top,nprobes)

	gene.selection <- match.arg(gene.selection,c("pairwise","common"))

	cn <- colnames(x)
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))
	topindex <- nprobes-top+1
	if(gene.selection=="pairwise") {
#		Distance measure is mean of top squared deviations for each pair of arrays
		for (i in 2:(nsamples))
		for (j in 1:(i-1))
			dd[i,j]=sqrt(mean(sort.int((x[,i]-x[,j])^2,partial=topindex)[topindex:nprobes]))
	} else {
#		Same genes used for all comparisons
		s <- rowMeans((x-rowMeans(x))^2)
		q <- quantile(s,p=(topindex-1.5)/(nprobes-1))
		x <- x[s>=q,]
		for (i in 2:(nsamples))
			dd[i,1:(i-1)]=sqrt(colMeans((x[,i]-x[,1:(i-1),drop=FALSE])^2))
	}

	a1 <- cmdscale(as.dist(dd),k=ndim)
	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top,gene.selection=gene.selection))
	mds$x <- a1[,dim.plot[1]]
	mds$y <- a1[,dim.plot[2]]
	plotMDS(mds,labels=labels,col=col,cex=cex,xlab=xlab,ylab=ylab,...)
}
