normexp.fit.control <- function(x, status=NULL, negctrl="negative", regular="regular", robust=FALSE)
#	Estimate normexp parameters using negative control probes
#  Wei Shi and Gordon Smyth
#  Created 17 April 2009. Last modified 19 April 2010.
{
	if(is(x, "EListRaw")){
		if(!is.null(status))
			s <- status
		else
			if(is.null(s <- x$genes$Status))
				stop("Probe status can not be found!")
		x$E <- as.matrix(x$E)
		xr <- x$E[tolower(s)==tolower(regular),,drop=FALSE]
		xn <- x$E[tolower(s)==tolower(negctrl),,drop=FALSE]
	}
	else {
		if(is.null(status)) stop("Please provide probe status!")
		x <- as.matrix(x)
		xr <- x[tolower(status)==tolower(regular),,drop=FALSE]
		xn <- x[tolower(status)==tolower(negctrl),,drop=FALSE]
	}
	
	if(robust) {
		require(MASS)
		narrays <- ncol(xn)
		m <- s <- rep(0,narrays)
	#	Robustness is judged on the log-scale, assumed normal
		for (j in 1:ncol(xn)) {
			h <- huber(log(xn[,j]))
			m[j] <- h$mu
			s[j] <- h$s
		}
	#	Means and standard deviation are converted back to log-normal
		mu <- exp(m+s^2/2)
		omega <- exp(s^2)
		sigma <- sqrt(omega*(omega-1))*exp(m)
	} else {
		mu <- colMeans(xn,na.rm=TRUE)
		sigma <- sqrt(rowSums((t(xn)-mu)^2,na.rm=TRUE)/(nrow(xn)-1))
	}
	alpha <- pmax(colMeans(xr,na.rm=TRUE)-mu,10)
	cbind(mu=mu,logsigma=log(sigma),logalpha=log(alpha))
}

nec <- function(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE)
#	Normexp background correction aided by negative controls.
#	Wei Shi
#	Created 27 September 2010.
{
	if(is(x,"EListRaw") && !is.null(x$Eb)) x$E <- x$E-x$Eb
	normexp.par <- normexp.fit.control(x, status=status, negctrl=negctrl, regular=regular, robust=robust)
	if(is(x, "EListRaw")) {
		for(i in 1:ncol(x))
			x$E[, i] <- normexp.signal(normexp.par[i, ], x$E[, i])
		x$E <- x$E + offset
	} else {
		x <- as.matrix(x)
		for(i in 1:ncol(x))
		x[, i] <- normexp.signal(normexp.par[i, ], x[, i])
		x <- x + offset
	}
	x
}

neqc <- function(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE, ...)
#	Normexp background correction and quantile normalization using control probes
#	Wei Shi and Gordon Smyth
#	Created 17 April 2009. Last modified 27 September 2010.
{
	x.bg <- nec(x,status,negctrl,regular,offset,robust)
	if(is(x.bg, "EListRaw")) {
		y <- normalizeBetweenArrays(x.bg, method="quantile", ...)
		if(is.null(status))
			status <- y$genes$Status
		y <- y[tolower(status) == tolower(regular), ]
		y$genes$Status <- NULL
	} else {
		x.bg <- as.matrix(x.bg)
		y <- log2(normalizeBetweenArrays(x.bg, method="quantile", ...))
		y <- y[tolower(status) == tolower(regular), ]
	}
	y
}
