#  M-A PLOTS

plotMA <- function(MA, array=1, xlab="A", ylab="M", main=colnames(MA$M)[array], xlim, ylim, status, values, pch, col, cex, legend=TRUE, ...) {
#	MA-plot with color coding for controls
#	Gordon Smyth  7 April 2003.
#	Revised by James Wettenhall  27 June 2003.
#	Last modified by GKS  10 Nov 2003.

	if(is(MA,"RGList")) {
		MA <- MA.RG(MA[,array])
		array <- 1
	}
	x <- as.matrix(MA$A)[,array]
	y <- as.matrix(MA$M)[,array]
	if(is.null(x) || is.null(y)) stop("No data to plot")
	if(missing(xlim)) xlim <- range(x,na.rm=TRUE)
	if(missing(ylim)) ylim <- range(y,na.rm=TRUE)
	if(missing(status)) status <- MA$genes$Status
	plot(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,type="n",...)
	if(is.null(status) || all(is.na(status))) {
		if(missing(pch)) pch=16
		if(missing(cex)) cex=0.2
		points(x,y,pch=pch[[1]],cex=cex[1])
	} else {
		if(missing(values)) {
			if(is.null(attr(status,"values")))
				values <- names(sort(table(status),decreasing=TRUE))
			else
				values <- attr(status,"values")
		}
#		Non-highlighted points
		sel <- !(status %in% values)
		nonhi <- any(sel)
		if(nonhi) points(x[sel],y[sel],pch=16,cex=0.2)

		nvalues <- length(values)
		if(missing(pch)) {
			if(is.null(attr(status,"pch")))
				pch <- rep(16,nvalues)
			else
				pch <- attr(status,"pch")
		}
		if(missing(cex)) {
			if(is.null(attr(status,"cex"))) {
				cex <- rep(1,nvalues)
				if(!nonhi) cex[1] <- 0.2
			} else
				cex <- attr(status,"cex")
		}
		if(missing(col)) {
			if(is.null(attr(status,"col"))) {
				col <- nonhi + 1:nvalues
			} else
				col <- attr(status,"col")
		}
		pch <- rep(pch,length=nvalues)
		col <- rep(col,length=nvalues)
		cex <- rep(cex,length=nvalues)
		for (i in 1:nvalues) {
			sel <- status==values[i]
			points(x[sel],y[sel],pch=pch[[i]],cex=cex[i],col=col[i])
		}
		if(legend) {
			if(is.list(pch))
				legend(x=xlim[1],y=ylim[2],legend=values,fill=col,col=col,cex=0.9)
			else
				legend(x=xlim[1],y=ylim[2],legend=values,pch=pch,,col=col,cex=0.9)
		}
	}
	invisible()
}

