plotWithHighlights <- function(x, y, status=NULL, values=NULL, pch=16, col=NULL, cex=1, legend="topleft", pch.bg=16, col.bg="black", cex.bg=0.3, ...)
#	Scatterplot with color coding for special points

#	Replaces the earlier function .plotMAxy, which in turn was based the original plotMA
#	created by Gordon Smyth 7 April 2003 and modified by James Wettenhall 27 June 2003.

#	Gordon Smyth
#	Last modified 17 April 2014.
{
#	If no status information, just plot all points normally
	if(is.null(status) || all(is.na(status))) {
		plot(x,y,pch=pch.bg,col=col.bg,cex=cex.bg,...)
		return(invisible())
	}
#	From here, status is not NULL and not all NA

#	Check values
	if(is.null(values)) {

#		Check for values and graphics parameters set as attributes by controlStatus()
		if(!is.null(attr(status,"values"))) {
			values <- attr(status,"values")
			if(!is.null(attr(status,"pch"))) pch <- attr(status,"pch")
			if(!is.null(attr(status,"col"))) col <- attr(status,"col")
			if(!is.null(attr(status,"cex"))) cex <- attr(status,"cex")
		}

#		Default is to set the most frequent status value as background, and to highlight other status values in decreasing order of frequency
		if(is.null(values)) {
			status.values <- names(sort(table(status),decreasing=TRUE))
			status <- as.character(status)
			values <- status.values[-1]
		}
	}

#	If no values, then just plot all points normally
	nvalues <- length(values)
	if(nvalues==0L) {
		plot(x,y,pch=pch.bg,col=col.bg,cex=cex.bg,...)
		return(invisible())
	}
#	From here, values has positive length

#	Setup plot axes
	plot(x,y,type="n",...)

#	Plot background (non-highlighted) points
	bg <- !(status %in% values)
	bg[is.na(bg)] <- TRUE
	nonhi <- any(bg)
	if(nonhi) points(x[bg],y[bg],pch=pch.bg[1],col=col.bg[1],cex=cex.bg[1])

#	Check graphical parameters for highlighted points
	pch <- rep_len(unlist(pch),length.out=nvalues)
	cex <- rep_len(unlist(cex),length.out=nvalues)
	if(is.null(col)) col <- nonhi + 1L:nvalues
	col <- rep_len(unlist(col),length.out=nvalues)

#	Plot highlighted points
	for (i in 1:nvalues) {
		sel <- status==values[i]
		points(x[sel],y[sel],pch=pch[i],cex=cex[i],col=col[i])
	}

#	Check legend
	if(is.logical(legend)) {
		legend.position <- "topleft"
	} else {
		legend.position <- as.character(legend)
		legend <- TRUE
	}
	legend.position <- match.arg(legend.position,c("bottomright","bottom","bottomleft","left","topleft","top","topright","right","center"))

#	Plot legend
	if(legend) {
		if(nonhi) {
#			Include background value in legend
			bg.value <- unique(status[bg])
			if(length(bg.value) > 1) bg.value <- "Other"
			values <- c(bg.value,values)
			pch <- c(pch.bg,pch)
			col <- c(col.bg,col)
			cex <- c(cex.bg,cex)
		}
		h <- cex>0.5
		cex[h] <- 0.5+0.8*(cex[h]-0.5)
		legend(legend.position,legend=values,pch=pch,,col=col,cex=0.9,pt.cex=cex)
	}

	invisible()
}
