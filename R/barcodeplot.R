##  BARCODEPLOT.R

barcodeplot <- function (statistics, index, index2 = NULL, labels = c("Up", "Down"), quantiles = c(-1, 1), col.bars = NULL, worm.window.width = 0, ...)
#	Barcode plot of one or two gene sets.
#	Gordon Smyth, Di Wu and Yifang Hu
#  20 October 2008.  Last revised 20 Jan 2014.
{
	# calculate sliding window sum
	slide.win.sum <- function(y, width) {
		n <- length(y)
		if(width > n) worm <- NULL else if(width <= n & width > 0) {
			cs <- c(0, cumsum(y))
			worm <- cs[(width+1):(n+1)] - cs[1:(n-width+1)]
			worm <- c(rep(0, ceiling(width/2-1)), worm, rep(0, floor(width/2)))
		}
		worm
	}
	
	# rescale x to new range
	rescale <- function(x, newrange) {
		r <- range(x)
		newrange[1] + (x-r[1]) / (r[2]-r[1]) * (newrange[2] - newrange[1])
	}

	TWO <- !is.null(index2)
	rankstat <- rank(statistics, na.last = TRUE)
	r <- rankstat[index]
	if (TWO) r2 <- rankstat[index2]
	isna <- is.na(statistics)
	n <- sum(!isna)
	
	if (any(isna)) {
		rankstat <- rankstat[1:n]
		r <- r[r <= n]
		if (TWO) r2 <- r2[r2 <= n]
	}
	
	if (is.null(col.bars)) 
		if (is.null(index2)) 
			col.bars = c("black", "black")
		else
			col.bars = c("red", "blue")
	quantiles <- sort(quantiles)

	# worm plot setting
	ylim.worm <- ylim <- c(-1, 1)
	ylab.worm <- xlab.worm <- ""
	if(!TWO) ylim.worm <- c(0, 1)
	
	if(worm.window.width == round(worm.window.width) & worm.window.width > 0) {
		ylim.worm <- c(-2.1, 2.1)
		if(!TWO) ylim.worm <- c(0, 2.1)
		xlab.worm <- "statistics"
	}
	
	ylim[2] <- ylim[2] + 0.5
	if (!is.null(index2)) ylim[1] <- ylim[1] - 0.5
	
	n <- length(statistics)
	
	plot(1:n, xlim = c(0, n), ylim = ylim.worm, type = "n", axes = FALSE, xlab = xlab.worm, ylab = ylab.worm, ...)
	
	npos <- sum(statistics > quantiles[2])
	nneg <- sum(statistics < quantiles[1])
	
	rect.yb <- -0.5
	if(!TWO) rect.yb <- 0
	
	rect(npos + 0.5, rect.yb, n - nneg + 0.5, 0.5, col = "lightgray", border = NA)
	if (npos) rect(0.5, rect.yb, npos + 0.5, 0.5, col = "pink", border = NA)
	if (nneg) rect(n - nneg + 0.5, rect.yb, n + 0.5, 0.5, col = "lightblue", border = NA)

	r <- n + 1 - r
	lwd <- 50/length(r)
	lwd <- min(2, lwd)
	lwd <- max(0.1, lwd)
	
	barlim <- ylim[2] - c(1.5, 0.5)  
	segments(r, barlim[1], r, barlim[2], lwd = lwd, col = col.bars[1])

	if(TWO) {
		r2 <- n + 1 - r2
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
	
	# create worm
	if(worm.window.width == round(worm.window.width) & worm.window.width > 0) {
		set <- matrix(0, nrow = length(statistics), ncol = 2)
		set[index, 1] <- 1    
		if(TWO) set[index2, 2] <- -1
		set <- set[order(statistics, decreasing = T),]
		
		# calculate enrichment
		ave.enrich1 <- length(index)/n
		worm1 <- (slide.win.sum(set[,1], worm.window.width)/worm.window.width)/ave.enrich1
		
		if(TWO) {
			ave.enrich2 <- length(index2)/n
			worm2 <- (slide.win.sum(set[,2], worm.window.width)/worm.window.width)/ave.enrich2
		}
		
		# rescale worm
		max.worm1 <- max(abs(worm1))
		worm1.scale <- rescale(worm1, c(1.1, 2.1))
		
		if(TWO) {
			max.worm2 <- max(abs(worm2))
			worm2.scale <- rescale(worm2, c(-2.1, -1.1))
		}
		
		# label statistics on x axis
		at.prob <- (0:10)/10 
		prob <- (10:0)/10
		axis(at = quantile(1:n, p = at.prob), side = 1, cex.axis = 0.8, las = 2, labels = format(quantile(statistics, p = prob), digits = 1))

		# plot worms
		worm.zero <- c(1:ceiling(worm.window.width/2-1), n:ceiling(n-worm.window.width/2+1))

		if(!TWO) {
			lines(x = c(1:n)[-worm.zero], y = worm1.scale[-worm.zero], col = col.bars[1], lwd = 2)
			axis(side = 2, at = c(1.1, 2.1), cex.axis = 0.8, labels = c(0, format(max.worm1, digits = 1)))
			axis(side = 2, labels = "Enrichment", at = 1.6, padj = -0.6, tick = FALSE, cex.axis = 0.8)
		}
	
		if(TWO) {
			lines(x = c(1:n)[-worm.zero], y = worm1.scale[-worm.zero], col = col.bars[1], lwd = 2)
			lines(x = c(1:n)[-worm.zero], y = worm2.scale[-worm.zero], col = col.bars[2], lwd = 2)
			axis(side = 2, at = c(1.1, 2.1), cex.axis = 0.7, labels = c(0, format(max.worm1, digits = 1)))
			axis(side = 2, at = c(-1.1, -2.1), cex.axis = 0.7, labels =  c(0, format(max.worm2, digits = 1)))
			axis(side = 2, labels = "Enrichment", at = -1.6, tick = FALSE, padj = -0.6, cex.axis = 0.7)
			axis(side = 2, labels = "Enrichment", at = 1.6, tick = FALSE, padj = -0.6, cex.axis = 0.7)
		}
	}

	invisible()
}