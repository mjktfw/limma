#  M-A PLOTS

plotMA <- function(MA, array=1, pch=16, status=NULL,
	values=c("gene","blank","buffer","utility","negative","calibration","ratio"),
	col=c("black","yellow","orange","pink","brown","blue","red"),
	cex=c(0.1,0.6,0.6,0.6,0.6,0.6,0.6)) {
#	MA-plot with color coding for controls
#	Gordon Smyth  7 April 2003.  Last modified 27 June 2003.
#	Revised by James Wettenhall  27 June 2003.

	x <- MA$A[,array]
	y <- MA$M[,array]
	if (is.list(pch)) 
		isListPCH <- TRUE 
	else 
		isListPCH <- FALSE
	pch <- as.list(pch)
	plot(x,y,xlab="A",ylab="M",main=colnames(MA$M)[array],type="n")
	if(is.null(status))
		points(x,y,pch=pch[[1]],cex=cex[1])
	else {
		nvalues <- length(values)
		if (length(pch) < nvalues)
			pch <- rep(pch,length=nvalues)
		col <- rep(col,nvalues)
		cex <- rep(cex,nvalues)
		for (i in 1:nvalues) {
			sel <- status==values[i]
			points(x[sel],y[sel],pch=pch[[i]],cex=cex[i],col=col[i])
		}
		if (isListPCH)
			legend(min(x,na.rm=TRUE),fill=col,max(y,na.rm=TRUE),legend=values,col=col,cex=0.9)
		else
			legend(min(x,na.rm=TRUE),pch=unlist(pch),max(y,na.rm=TRUE),legend=values,col=col,cex=0.9)
	}
	invisible()
}

