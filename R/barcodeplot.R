##  BARCODEPLOT.R

barcodeplot <- function (statistics, index, index2 = NULL, labels = c("Up", "Down"), quantiles = c(-1, 1), col.bars = NULL, worm=TRUE, span.worm=0.45, ...)
#	Barcode plot of one or two gene sets.
#	Gordon Smyth, Di Wu and Yifang Hu
#	20 October 2008.  Last revised 4 Feb 2014.
{
#	Are there up and down sets?
	TWO <- !is.null(index2)

#	Convert indexes to logical
	set1 <- set2 <- rep.int(FALSE,length(statistics))
	names(set1) <- names(set2) <- names(statistics)
	set1[index] <- TRUE
	if(TWO) set2[index2] <- TRUE

#	Sort statistics and indexes
	ostat <- order(statistics, na.last = TRUE, decreasing=TRUE)
	statistics <- statistics[ostat]
	set1 <- set1[ostat]
	if(TWO) set2 <- set2[ostat]

#	Trim off missing values
	n <- sum(!is.na(statistics))
	if(n==0L) {
		message("No valid statistics")
		return(invisible())
	}
	if (n < length(statistics)) {
		statistics <- statistics[1:n]
		set1 <- set1[1:n]
		if (TWO) set2 <- set2[1:n]
	}

#	Convert indexes to integer
	r <- which(set1)
	if (TWO) {
		r2 <- which(set2)
		if(!length(r2)) TWO <- FALSE
	}

#	Check there is something to plot
	if(!length(r))
		if(TWO) {
			r <- r2
			set1 <- set2
			TWO <- FALSE
			message("Using index2 as primary set")
		} else {
			message("No selected genes to plot")
			return(invisible())
		}

#	Check other arguments
	quantiles <- sort(quantiles)
	if (is.null(col.bars)) 
		if (TWO) 
			col.bars = c("red", "blue")
		else
			col.bars = c("black", "black")

	# worm plot setting
	ylim.worm <- ylim <- c(-1, 1)
	ylab.worm <- ""
	xlab.worm <- "statistics"
	
	if(!TWO) ylim.worm <- c(0, 1)
	
	if(worm) {
		ylim.worm <- c(-2.1, 2.1)
		if(!TWO) ylim.worm <- c(0, 2.1)
	}
	
	ylim[2] <- ylim[2] + 0.5
	if (TWO) ylim[1] <- ylim[1] - 0.5
	
	plot(1:n, xlim = c(0, n), ylim = ylim.worm, type = "n", axes = FALSE, xlab = xlab.worm, ylab = ylab.worm, ...)
	
	npos <- sum(statistics > quantiles[2])
	nneg <- sum(statistics < quantiles[1])

	rect.yb <- -0.5
	if(!TWO) rect.yb <- 0
	
	rect(npos + 0.5, rect.yb, n - nneg + 0.5, 0.5, col = "lightgray", border = NA)
	if (npos) rect(0.5, rect.yb, npos + 0.5, 0.5, col = "pink", border = NA)
	if (nneg) rect(n - nneg + 0.5, rect.yb, n + 0.5, 0.5, col = "lightblue", border = NA)

	lwd <- 50/length(r)
	lwd <- min(2, lwd)
	lwd <- max(0.1, lwd)

	barlim <- ylim[2] - c(1.5, 0.5)  
	segments(r, barlim[1], r, barlim[2], lwd = lwd, col = col.bars[1])

	if(TWO) {
		lwd2 <- 50/length(r2)
		lwd2 <- min(2, lwd2)
		lwd2 <- max(0.1, lwd2)
		barlim2 <- ylim[1] + c(0.5, 1.5)
		segments(r2, barlim2[1], r2, barlim2[2], lwd = lwd2, col = col.bars[2])
	}
	
	lab.at <- 0
	if(!TWO) lab.at <- 0.5
	axis(side = 2, at = lab.at, padj = 3, cex.axis = 0.85, labels = labels[1], tick = FALSE)
	axis(side = 4, at = lab.at, padj = -3, cex.axis = 0.85, labels = labels[2], tick = FALSE)
	
	# label statistics on x axis
	prob <- (10:0)/10
	axis(at = seq(1,n,len=11), side = 1, cex.axis = 0.7, las = 2, labels = format(quantile(statistics, p = prob), digits = 1))

	# create worm
	if(worm) {
		# rescale x to new range
		rescale <- function(x, newrange, oldrange=range(x)) {
			newrange[1] + (x-oldrange[1]) / (oldrange[2]-oldrange[1]) * (newrange[2] - newrange[1])
		}

		# calculate enrichment
		ave.enrich1 <- length(r)/n
		worm1 <- tricubeMovingAverage(set1,span=span.worm)/ave.enrich1
	
		if(TWO) {
			ave.enrich2 <- length(r2)/n
			worm2 <- tricubeMovingAverage(-set2,span=span.worm)/ave.enrich2
		}
	
		# rescale worm
		max.worm1 <- max(worm1)
		r.worm1 <- c(0,max.worm1)
		worm1.scale <- rescale(worm1, newrange=c(1.1,2.1), oldrange=r.worm1)
	
		if(TWO) {
			min.worm2 <- min(worm2)
			r.worm2 <- c(min.worm2,0)
			worm2.scale <- rescale(worm2, newrange=c(-2.1,-1.1), oldrange=r.worm2)
		}

		# plot worms

		if(!TWO) {
			lines(x = 1:n, y = worm1.scale, col = col.bars[1], lwd = 2)
			abline(h = rescale(1,newrange=c(1.1,2.1),oldrange=r.worm1), lty=2)
			axis(side = 2, at = c(1.1, 2.1), cex.axis = 0.8, labels = c(0, format(max.worm1, digits = 2)))
			axis(side = 2, labels = "Enrichment", at = 1.6, padj = -0.6, tick = FALSE, cex.axis = 0.8)
		}
	
		if(TWO) {
			lines(x = 1:n, y = worm1.scale, col = col.bars[1], lwd = 2)
			abline(h = rescale(1,newrange=c(1.1,2.1),oldrange=r.worm1), lty=2)
			lines(x = 1:n, y = worm2.scale, col = col.bars[2], lwd = 2)
			abline(h = rescale(-1,newrange=c(-2.1,-1.1),oldrange=r.worm2), lty=2)
			axis(side = 2, at = c(1.1, 2.1), cex.axis = 0.7, labels = c(0, format(max.worm1, digits = 2)))
			axis(side = 2, at = c(-1.1, -2.1), cex.axis = 0.7, labels =  c(0, format(-min.worm2, digits = 2)))
			axis(side = 2, labels = "Enrichment", at = -1.6, tick = FALSE, padj = -0.6, cex.axis = 0.7)
			axis(side = 2, labels = "Enrichment", at = 1.6, tick = FALSE, padj = -0.6, cex.axis = 0.7)
		}
	}

	invisible()
}
