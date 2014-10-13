##  PLOTMDS.R

#	Class to hold multidimensional scaling output
setClass("MDS",representation("list"))

setMethod("show","MDS",function(object) {
	cat("An object of class MDS\n")
	print(unclass(object))
})

plotMDS <- function(x,...) UseMethod("plotMDS")

plotMDS.MDS <- function(x,labels=NULL,pch=NULL,cex=1,dim.plot=x$dim.plot,xlab=paste("Dimension",dim.plot[1]),ylab=paste("Dimension",dim.plot[2]),...)
#	Method for MDS objects
#	Create a new plot using MDS coordinates or distances previously created
#	Gordon Smyth and Yifang Hu
#	21 May 2011.  Last modified 26 June 2014
{
#	Check labels
	if(is.null(labels) & is.null(pch)) {
		labels <- colnames(x$distance.matrix)
		if(is.null(labels)) labels <- 1:length(x$x)
	}

#	Are new dimensions requested?
	if(!all(dim.plot==x$dim.plot)) {
		ndim <- max(dim.plot)
		if(ndim > ncol(x$cmdscale.out)) x$cmdscale.out <- cmdscale(as.dist(x$distance.matrix),k=ndim)
		x$x <- x$cmdscale.out[,dim.plot[1]]
		x$y <- x$cmdscale.out[,dim.plot[2]]
	}

#	Make the plot
	if(is.null(labels)){
#		Plot symbols instead of text
		plot(x$x, x$y, pch = pch, xlab = xlab, ylab = ylab, cex = cex, ...)
	} else {
#		Plot text.  Need to estimate width of labels in plot coordinates.
#		Estimate will be ok for default plot width, but maybe too small for smaller plots.
		labels <- as.character(labels)
		StringRadius <- 0.01*cex*nchar(labels)
		left.x <- x$x-StringRadius
		right.x <- x$x+StringRadius
		plot(c(left.x, right.x), c(x$y, x$y), type = "n", xlab = xlab, ylab = ylab, ...)
		text(x$x, x$y, labels = labels, cex = cex, ...)
	}

	return(invisible(x))
}

plotMDS.default <- function(x,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),ndim=max(dim.plot),gene.selection="pairwise",xlab=paste("Dimension",dim.plot[1]),ylab=paste("Dimension",dim.plot[2]),...)
#	Multi-dimensional scaling with top-distance
#	Di Wu and Gordon Smyth
#	19 March 2009.  Last modified 26 June 2014
{
#	Check x
	x <- as.matrix(x)
	nsamples <- ncol(x)
	if(nsamples < 3) stop("Need at least 3 columns")
	cn <- colnames(x)
#	Remove rows with missing or Inf values
	bad <- rowSums(is.finite(x)) < nsamples
	if(any(bad)) x <- x[!bad,,drop=FALSE]
	nprobes <- nrow(x)

#	Check top
	top <- min(top,nprobes)

#	Check labels and pch
	if(is.null(pch) & is.null(labels)) {
		labels <- colnames(x)
		if(is.null(labels)) labels <- 1:nsamples
	}
	if(!is.null(labels)) labels <- as.character(labels)

#	Check dim
	if(ndim < 2) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("Too few samples")
	if(nprobes < ndim) stop("Too few rows")

#	Check gene.selection
	gene.selection <- match.arg(gene.selection,c("pairwise","common"))

#	Distance matrix from pairwise leading fold changes
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

#	Multi-dimensional scaling
	a1 <- cmdscale(as.dist(dd),k=ndim)

#	Make MDS object and call plotMDS method
	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top,gene.selection=gene.selection))
	mds$x <- a1[,dim.plot[1]]
	mds$y <- a1[,dim.plot[2]]
	plotMDS(mds,labels=labels,pch=pch,cex=cex,xlab=xlab,ylab=ylab,...)
}
