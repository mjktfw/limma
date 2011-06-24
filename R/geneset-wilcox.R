##  GENESET.R

geneSetTest <- function(selected,statistics,alternative="mixed",type="auto",ranks.only=TRUE,nsim=10000)
#	Gene set test using either Wilcox test or simulation.
#	Gordon Smyth
#	3 September 2004. Last modified 21 July 2006.
{
	alternative <- match.arg(alternative,c("mixed","either","down","up","less","greater"))
	if(alternative=="two.sided") alternative <- "either"
	if(alternative=="less") alternative <- "down"
	if(alternative=="greater") alternative <- "up"
	type <- match.arg(tolower(type),c("auto","t","f"))
	allsamesign <- all(statistics >= 0) || all(statistics <= 0)
	if(type=="auto") {
		if(allsamesign)
			type <- "f"
		else
			type <- "t"
	}
	if(type=="f" & alternative!="mixed") stop("Only alternative=\"mixed\" is possible with F-like statistics.")
	if(alternative=="mixed") statistics <- abs(statistics)
	if(alternative=="down") {
		statistics <- -statistics
		alternative <- "up"
	}
	if(ranks.only) {
#		The test statistic is the mean rank of the selected statistics
#		and the p-value is obtained explicitly from the Wilcox test
		if(alternative=="either")
			wilc.alt <- "two.sided"
		else
			wilc.alt <- "greater"
		x <- y <- NULL
		if(is.logical(selected)) {
			x <- statistics[selected]
			y <- statistics[!selected]
		}
		if(is.numeric(selected)) {
			x <- statistics[selected]
			y <- statistics[-selected]
		}
		if(is.character(selected)) {
			nam <- names(statistics)
			if(is.null(nam)) stop("selected is character but elements of statistics are not named")
			selected <- is.element(nam,selected)
			x <- statistics[selected]
			y <- statistics[!selected]
		}
		return(wilcox.test(x,y,alternative=wilc.alt,conf.int=FALSE)$p.value)
	} else {
#		The test statistic is the mean of the selected statistics
#		and the p-value is obtained by random permutation
		ssel <- statistics[selected]
		ssel <- ssel[!is.na(ssel)]
		nsel <- length(ssel)
		if(nsel==0) return(1)
		stat <- statistics[!is.na(statistics)]
		msel <- mean(ssel)
		if(alternative=="either")
			posstat <- abs
		else
			posstat <- function(x) x
		msel <- posstat(msel)
		ntail <- 0
		for (i in 1:nsim) if(posstat(mean(sample(stat,nsel))) >= msel) ntail <- ntail+1
		return(ntail/nsim)
	}
}

wilcoxGST <- function(selected,statistics,alternative="mixed")
#	Mean-rank gene set test using Wilcox test.
#	Gordon Smyth
#	27 July 2009.  Last modified 28 July 2009.
geneSetTest(selected=selected,statistics=statistics,alternative=alternative,type="t",ranks.only=TRUE)

barcodeplot <- function(selected,statistics,labels=c("Up","Down"),...)
#	Barcode plot for gene set test
#	Gordon Smyth and Di Wu
#  20 October 2008. Last revised 30 Sep 2009.
{
	isna <- is.na(statistics)
	if(any(isna)) {
		if(length(selected)==length(statistics)) {
			selected <- selected[!isna]
		} else {
			selected <- as.integer(selected)
			selected <- selected[!isna[selected]]
		}
		statistics <- statistics[!isna]
	}

	n <- length(statistics)
	plot(1:n,xlim=c(0,n),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="",...)
	npos <- sum(statistics > 1)
	nneg <- sum(statistics < -1)
	rect(npos+0.5,0,n-nneg+0.5,1,col="lightgray",border=NA)
	if(npos) rect(0.5,0,npos+0.5,1,col="pink",border=NA)
	if(nneg) rect(n-nneg+0.5,0,n+0.5,1,col="lightgreen",border=NA)
	r <- n+1-rank(statistics)[selected]
	lwd <- 50/length(r)
	lwd <- min(2,lwd)
	lwd <- max(0.1,lwd)
	segments(r,0,r,1,lwd=lwd)
#	rect(0.5,0,n+0.5,1,border="blue")
	mtext(labels[1],side=2,line=-1,col="gray")
	mtext(labels[2],side=4,line=-1,col="gray")
	invisible()
}

barcodeplot2 <- function (selected,statistics,selected2=NULL,labels=c("Up","Down"), ...) 
#	Barcode plot for gene set test. 
#  Plots both up and down sets in one figure.
#	Gordon Smyth and Di Wu
#  20 October 2008. Last revised 30 Sep 2009.
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
	n <- length(statistics)

	plot(1:n,xlim=c(0,n),ylim=c(-1,1),type="n",axes=FALSE,xlab="",ylab="",...)
	npos <- sum(statistics > 1)
	nneg <- sum(statistics < -1)
	rect(npos+0.5,-0.5,n-nneg+0.5,0.5,col="lightgray",border=NA)
	if(npos) rect(0.5,-0.5,npos+0.5,0.5,col="pink",border=NA)
	if(nneg) rect(n-nneg+0.5,-0.5,n+0.5,0.5,col="lightblue",border=NA)

	rankstat <- rank(statistics)
	r <- n+1-rankstat[selected]
	lwd <- 50/length(r)
	lwd <- min(2,lwd)
	lwd <- max(0.1,lwd)
	segments(r,0,r,1,lwd=lwd,col="red")

	if(!is.null(selected2)) {
		r2 <- n+1-rankstat[selected2]
		lwd2 <- 50/length(r2)
		lwd2 <- min(2,lwd2)
		lwd2 <- max(0.1,lwd2)
		segments(r2,-1,r2,0,lwd=lwd2,col="blue")
	}

	mtext(labels[1],side=2,line=-0.5,col="black")
	mtext(labels[2],side=4,line=-1,col="black")
	invisible()
}

