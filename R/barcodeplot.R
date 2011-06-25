##  BARCODEPLOT.R

barcodeplot <- function(statistics,selected,selected2=NULL,labels=c("Up","Down"),quantiles=c(-1,1),col.bars=NULL,offset.bars=!is.null(selected2), ...) 
#	Barcode plot of one or two gene sets. 
#	Gordon Smyth and Di Wu
#  20 October 2008.  Last revised 25 June 2011.
{
	isna <- is.na(statistics)
	if (any(isna)) {
		if (length(selected) == length(statistics)) {
			selected <- selected[!isna]
		}
		else {
			selected <- as.integer(selected)
			selected <- selected[!isna[selected]]
		}
		if(!is.null(selected2))
			if (length(selected2) == length(statistics)) {
				selected2 <- selected2[!isna]
			}
			else {
				selected2 <- as.integer(selected2)
				selected2 <- selected[!isna[selected2]]
			}
		statistics <- statistics[!isna]
	}
	if(is.null(col.bars))
		if(is.null(selected2))
			col.bars=c("black","black")
		else
			col.bars=c("red","blue")
	quantiles <- sort(quantiles)

	ylim <- c(-1,1)
	if(offset.bars) {
		ylim[2] <- ylim[2]+0.5
		if(!is.null(selected2)) ylim[1] <- ylim[1]-0.5
	}
	n <- length(statistics)
	plot(1:n,xlim=c(0,n),ylim=ylim,type="n",axes=FALSE,xlab="",ylab="",...)
	npos <- sum(statistics > quantiles[2])
	nneg <- sum(statistics < quantiles[1])
	rect(npos+0.5,-0.5,n-nneg+0.5,0.5,col="lightgray",border=NA)
	if(npos) rect(0.5,-0.5,npos+0.5,0.5,col="pink",border=NA)
	if(nneg) rect(n-nneg+0.5,-0.5,n+0.5,0.5,col="lightblue",border=NA)

	rankstat <- rank(statistics)
	r <- n+1-rankstat[selected]
	lwd <- 50/length(r)
	lwd <- min(2,lwd)
	lwd <- max(0.1,lwd)
	barlim <- ylim[2]-c(1.5,0.5)
	segments(r,barlim[1],r,barlim[2],lwd=lwd,col=col.bars[1])

	if(!is.null(selected2)) {
		r2 <- n+1-rankstat[selected2]
		lwd2 <- 50/length(r2)
		lwd2 <- min(2,lwd2)
		lwd2 <- max(0.1,lwd2)
		barlim2 <- ylim[1]+c(0.5,1.5)
		segments(r2,barlim2[1],r2,barlim2[2],lwd=lwd2,col=col.bars[2])
	}

	mtext(labels[1],side=2,line=-0.5,col="black")
	mtext(labels[2],side=4,line=-1,col="black")
	invisible()
}

