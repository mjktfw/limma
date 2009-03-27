##  PLOTMDS.R

plotMDS <- function(x,top=500,labels=colnames(x),col=NULL,cex=1,dim.plot=c(1,2),ndim=max(dim.plot),...)
#	Multi-dimensional scaling with top-distance
#	Di Wu and Gordon Smyth
#	19 March 2009.  Last modified 26 March 2009.
{
	x <- as.matrix(x)
	if(is.null(labels)) {labels<-1:dim(x)[2]}

	nprobes <- nrow(x)
	nsamples <- ncol(x)
	if(ndim < 2) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("Two few samples")

#	Distance measure is mean of top squared deviations for each pair of arrays
	dd <- matrix(0,nrow=nsamples,ncol=nsamples)
	topindex <- nprobes-top+1
	for (i in 2:(nsamples))
	for (j in 1:(i-1))
		dd[i,j]=sqrt(mean(sort.int((x[,i]-x[,j])^2,partial=topindex)[topindex:nprobes]))

	a1 <- cmdscale(as.dist(dd),k=ndim)
	plot(a1[,dim.plot[1]],a1[,dim.plot[2]],type="n",xlab=paste("Dimension",dim.plot[1]),ylab=paste("Dimension",dim.plot[2]),...)
	text(a1[,dim.plot[1]],a1[,dim.plot[2]],labels=labels,col=col,cex=cex)
	return(invisible(dd))
}
