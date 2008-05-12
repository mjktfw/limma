#  BACKGROUND.R

#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG, method="subtract", offset=0, printer=RG$printer, verbose=TRUE) {
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003.  Last modified 12 May 2008.

	if(!is.list(RG) && is.vector(RG)) RG <- as.matrix(RG)
	if(is.matrix(RG)) {
		method <- match.arg(method, c("none","normexp","saddle","neldermean","bfgs","rma","mcgee"))
		if(method=="normexp") method="saddle"
		if(method!="none") {
			for (j in 1:ncol(RG)) {
				x <- RG[,j]
				out <- normexp.fit.C(x,method=method)
				RG[,j] <- normexp.signal(out$par,x)
				if(verbose) cat("Corrected array",j,"\n")
			}
		}
#		rma={
#			require("affy")
#			RG <- apply(RG,2,bg.adjust)
#		})
		if(offset) RG <- RG+offset
		return(RG)
	}

	if(is.null(RG$Rb) != is.null(RG$Gb)) stop("Background values exist for one channel but not the other")
	method <- match.arg(method, c("none","subtract","half","minimum","movingmin","edwards","normexp","saddle","neldermean","bfgs","rma","mcgee"))
	if(is.null(RG$Rb) && is.null(RG$Gb)) method <- "none"
	switch(method,
	subtract={
		RG$R <- RG$R-RG$Rb
		RG$G <- RG$G-RG$Gb
	},
	half={
		RG$R <- pmax(RG$R-RG$Rb, 0.5)
		RG$G <- pmax(RG$G-RG$Gb, 0.5)
	},
	minimum={
		RG$R <- as.matrix(RG$R - RG$Rb)
		RG$G <- as.matrix(RG$G - RG$Gb)
		for (slide in 1:ncol(RG$R)) {
			i <- RG$R[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(RG$R[!i,slide],na.rm=TRUE)
				RG$R[i,slide] <- m/2
			}
			i <- RG$G[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(RG$G[!i,slide],na.rm=TRUE)
				RG$G[i,slide] <- m/2
			}
		}
	},
	movingmin={
		RG$R <- RG$R-ma3x3.spottedarray(RG$Rb,printer=printer,FUN=min,na.rm=TRUE)
		RG$G <- RG$G-ma3x3.spottedarray(RG$Gb,printer=printer,FUN=min,na.rm=TRUE)
	},
	edwards={
#		Log-linear interpolation for dull spots as in Edwards (2003).
#		The threshold values (delta) are chosen such that the number of
#		spots with (0 < R-Rb < delta) is f=10% of the number spots
#		with (R-Rb <= 0) for each channel and array.
#		Note slight change from Edwards (2003).
		one <- matrix(1,NROW(RG$R),1)
		delta.vec <- function(d, f=0.1) {
			quantile(d, mean(d<1e-16,na.rm=TRUE)*(1+f), na.rm=TRUE)
		}
		sub <- as.matrix(RG$R-RG$Rb)
		delta <- one %*% apply(sub, 2, delta.vec)
		RG$R <- ifelse(sub < delta, delta*exp(1-(RG$Rb+delta)/RG$R), sub)
		sub <- as.matrix(RG$G-RG$Gb)
		delta <- one %*% apply(sub, 2, delta.vec)
		RG$G <- ifelse(sub < delta, delta*exp(1-(RG$Gb+delta)/RG$G), sub)
	},
	normexp=,saddle=,neldermean=,bfgs=,rma=,mcgee={
		if(verbose) cat("Green channel\n")
		RG$G <- backgroundCorrect(RG$G-RG$Gb,method=method,verbose=verbose)
		if(verbose) cat("Red channel\n")
		RG$R <- backgroundCorrect(RG$R-RG$Rb,method=method,verbose=verbose)
	}
	)
	RG$Rb <- NULL
	RG$Gb <- NULL
	if(offset) {
		RG$R <- RG$R+offset
		RG$G <- RG$G+offset
	}
	new("RGList",unclass(RG))
}

ma3x3.matrix <- function(x,FUN=mean,na.rm=TRUE,...)
#	2-dimensional moving average for 3x3 blocks
#	Gordon Smyth
#	11 April 2004
{
#	Pad out x with NA so that original values have 8 neighbors
	d1 <- nrow(x)
	d2 <- ncol(x)
	y <- matrix(NA,d1+2,d2+2)
	y[1+(1:d1),1+(1:d2)] <- x

#	Index vector for original values
	i <- 1:length(y)
	dim(i) <- dim(y)
	i <- i[1+(1:d1),1+(1:d2)]
	dim(i) <- NULL

#	Rows are original obs, columns are neighbors
	x <- matrix(x,d1*d2,9)
	ry <- nrow(y)
	x[,1] <- y[i-ry-1]
	x[,2] <- y[i-ry]
	x[,3] <- y[i-ry+1]
	x[,4] <- y[i-1]
	x[,6] <- y[i+1]
	x[,7] <- y[i+ry-1]
	x[,8] <- y[i+ry]
	x[,9] <- y[i+ry+1]

	y <- apply(x,MARGIN=1,FUN=FUN,na.rm=na.rm,...)
	dim(y) <- c(d1,d2)
	y
}

ma3x3.spottedarray <- function(x,printer,FUN=mean,na.rm=TRUE,...)
#	Gordon Smyth
#	11 April 2004
{
	x <- as.matrix(x)
	narrays <- ncol(x)
	gr <- printer$ngrid.r
	gc <- printer$ngrid.c
	sr <- printer$nspot.r
	sc <- printer$nspot.c
	dim(x) <- c(sc, sr, gc, gr, narrays)
	x <- aperm(x, perm = c(2, 4, 1, 3, 5))
	dim(x) <- c(gr * sr, gc * sc, narrays)
	for (j in 1:narrays) x[,,j] <- ma3x3.matrix(x[,,j],FUN=FUN,na.rm=TRUE,...)
	dim(x) <- c(sr, gr, sc, gc, narrays)
	x <- aperm(x, perm = c(3, 1, 4, 2, 5))
	dim(x) <- c(sc*sr*gc*gr, narrays)
	x
}

