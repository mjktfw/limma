##  PLOTMDS.R

plotMDS <- function(x,top=500,labels=colnames(x),col=NULL,cex=1,dim.plot=c(1,2),ndim=max(dim.plot),gene.selection="pairwise",xlab=paste("Dimension",dim.plot[1]),ylab=paste("Dimension",dim.plot[2]),...)
#	Multi-dimensional scaling with top-distance
#	Di Wu and Gordon Smyth
#	19 March 2009.  Last modified 6 May 2011.
{
	x <- as.matrix(x)

#	Remove rows with missing or Inf values
	ok <- is.finite(x)
	if(!all(ok)) x <- x[apply(ok,1,all),]

	if(is.null(labels)) {labels<-1:dim(x)[2]}

	nprobes <- nrow(x)
	nsamples <- ncol(x)
	if(ndim < 2) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("Two few samples")

	gene.selection <- match.arg(gene.selection,c("pairwise","common"))

	dd <- matrix(0,nrow=nsamples,ncol=nsamples)
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
	xcord <- a1[,dim.plot[1]]
	ycord <- a1[,dim.plot[2]]
	plot(xcord,ycord,type="n",xlab=xlab,ylab=ylab,...)
	text(xcord,ycord,labels=labels,col=col,cex=cex)
	return(invisible(list(distance.matrix=dd,x=xcord,y=ycord)))
}

