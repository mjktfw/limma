#  WITHIN ARRAY NORMALIZATION

MA.RG <- function(object) {
#	Convert RGList to MAList
#	Gordon Smyth
#	2 March 2003.  Last revised 28 April 2003.

	R <- object$R
	G <- object$G

#	Background correction
	if(!is.null(object$Rb)) R <- R-object$Rb
	if(!is.null(object$Gb)) G <- G-object$Gb

#	Log
	R[R <= 0] <- NA
	G[G <= 0] <- NA
	R <- log(R,2)
	G <- log(G,2)
	
#	Minus and Add
	object$R <- object$G <- object$Rb <- object$Gb <- NULL
	object$M <- as.matrix(R-G)
	object$A <- as.matrix((R+G)/2)
	new("MAList",unclass(object))
}

normalizeWithinArrays <- function(object,layout,method="printtiploess",weights=object$weights,span=0.4,iterations=4,controlspots=NULL,df=5,robust="M") {
#	Sub-array loess normalization
#	Gordon Smyth
#	2 March 2003.  Last revised 15 June 2003.

	if(!is(object,"MAList")) object <- MA.RG(object)
	choices <- c("none","median","loess","printtiploess","composite","robustspline")
	method <- choices[pmatch(method,choices)]
	if(is.na(method)) warning("normalization method not recognized, defaulting to \"none\"")
	if(is.na(method) || method=="none") return(object)
	ngenes <- nrow(object$M)
	narrays <- ncol(object$M)
	switch(method,
		median = {
			for (j in 1:narrays) object$M[,j] <- object$M[,j] - median(object$M[,j],na.rm=TRUE)
		},
		loess = {
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				object$M[,j] <- residuals(loess(y~x,weights=w,span=0.75*span,na.action=na.exclude,degree=1,family="symmetric",trace.hat="approximate",iterations=iterations,surface="direct"))
			}
		},
		printtiploess = {
			for (j in 1:narrays) {
				ngr <- layout$ngrid.r
				ngc <- layout$ngrid.c
				nspots <- layout$nspot.r * layout$nspot.c
				spots <- 1:nspots
				for (gridr in 1:ngr)
				for (gridc in 1:ngc) {
					y <- object$M[spots,j]
					x <- object$A[spots,j]
					w <- weights[spots,j]
					object$M[spots,j] <- residuals(loess(y~x,weights=w,span=span,na.action=na.exclude,degree=1,family="symmetric",trace.hat="approximate",iterations=iterations,surface="direct"))
					spots <- spots + nspots
				}
			}
		},
		composite = {
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				fit <- loess(y~x,weights=w,span=span,subset=controlspots,na.action=na.exclude,degree=0,surface="direct",family="symmetric",trace.hat="approximate",iterations=iterations)
				global <- predict(fit,newdata=x)
				alpha <- (rank(x)-1) / sum(!is.na(x))
				ntips <- layout$ngrid.r * layout$ngrid.c
				nspots <- layout$nspot.r * layout$nspot.c
				spots <- 1:nspots
				for (tip in 1:ntips) {
					y <- object$M[spots,j]
					x <- object$A[spots,j]
					w <- weights[spots,j]
					local <- fitted(loess(y~x,weights=w,span=span,na.action=na.exclude,degree=1,family="symmetric",trace.hat="approximate",iterations=iterations,surface="direct"))
					object$M[spots,j] <- object$M[spots,j] - alpha[spots]*global[spots]-(1-alpha[spots])*local
					spots <- spots + nspots
				}
			}
		},
		robustspline = {
			for (j in 1:narrays)
				object$M[,j] <- normalizeRobustSpline(object$M[,j],object$A[,j],layout,df=df,method=robust)
		}
	)
	object
}

normalizeRobustSpline <- function(M,A,layout,df=5,method="M") {
#	Robust spline normalization
#	Gordon Smyth
#	27 April 2003.  Last revised 28 April 2003.

	require(MASS)
	require(splines)
	ngrids <- layout$ngrid.r * layout$ngrid.c
	nspots <- layout$nspot.r * layout$nspot.c
	col <- rainbow(ngrids,end=(ngrids-2)/ngrids)

#	Global splines
	O <- is.finite(M) & is.finite(A)
	X <- matrix(NA,ngrids*nspots,df)
	X[O,] <- ns(A[O],df=df,intercept=TRUE)
	x <- X[O,]
	y <- M[O]
	s <- summary(rlm(x,y,method=method))
	beta0 <- s$coefficients[,1]
	covbeta0 <- s$cov * s$stddev^2

#	Tip-wise splines
	beta <- array(0,c(ngrids,df))
	covbeta <- array(0,c(ngrids,df,df))
	spots <- 1:nspots
	for (i in 1:ngrids) {
		o <- O[spots]
		y <- M[spots][o]
		x <- X[spots,][o,]
		s <- summary(rlm(x,y,method=method))
		beta[i,] <- s$coefficients[,1]
		covbeta[i,,] <- s$cov * s$stddev^2
		spots <- spots + nspots
	}

#	Empirical Bayes estimates
	res.cov <- cov(beta) - apply(covbeta,c(2,3),mean)
	Sigma0 <- covbeta0 * max(0, sum(diag(res.cov)) / sum(diag(covbeta0)) )
#	Sigma0 <- covbeta0 * max(0,mean(eigen(solve(covbeta0,res.cov))$values))

#	Shrunk splines
	spots <- 1:nspots
	for (i in 1:ngrids) {
		beta[i,] <- beta0 + Sigma0 %*% solve(Sigma0+covbeta[i,,],beta[i,]-beta0)
		o <- O[spots]
		x <- X[spots,][o,]
		M[spots][o] <- M[spots][o] - x %*% beta[i,]
		M[spots][!o] <- NA
		spots <- spots + nspots
	}
	M
}


#  PRINTORDER

normalizeForPrintorder <- function(object,layout,start="topleft",method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order
#	Gordon Smyth
#	11 Mar 2002.  Last revised 18 June 2003.

	if(is.null(object$R) || is.null(object$G)) stop("List must contain components R and G")
	start <- match.arg(start,c("topleft","topright"))
	method <- match.arg(method,c("loess","plate"))
	po <- printorder(layout,start=start)$printorder
	nslides <- NCOL(object$R)
	for (i in 1:nslides) {
		RG <- normalizeForPrintorder.rg(R=object$R[,i],G=object$G[,i],printorder=po,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size)
		object$R[,i] <- RG$R
		object$G[,i] <- RG$G
	}
	object
}

normalizeForPrintorder.rg <- function(R,G,printorder,method="loess",separate.channels=FALSE,span=0.1,plate.size=32,plot=FALSE) {
#	Pre-normalize the foreground intensities for print order, given R and G for a single array.
#	Gordon Smyth
#	8 Mar 2002.  Last revised 18 June 2003.

	if(plot) ord <- order(printorder)
	Rf <- log(R,2)
	Gf <- log(G,2)
	Rf[is.infinite(Rf)] <- NA
	Gf[is.infinite(Gf)] <- NA
	if(!separate.channels) Af <- (Rf+Gf)/2
	method <- match.arg(method,c("loess","plate"))
	if(method=="plate") {
		# Correct for plate pack (usually four 384-well plates)
		plate <- 1 + (printorder-0.5) %/% plate.size
		hubermu <- function(...) huber(...)$mu
		if(separate.channels) {
			plate.mR <- tapply(Rf,plate,hubermu)
			plate.mG <- tapply(Gf,plate,hubermu)
			mR <- mG <- Rf
			for (i in 1:max(plate)) {
				mR[plate==i] <- plate.mR[i]
				mG[plate==i] <- plate.mG[i]
			}
			if(plot) {
				plot(printorder,Rf,xlab="Print Order",ylab="Log Intensity",type="n")
				points(printorder,Rf,pch=".",col="red")
				points(printorder,Gf,pch=".",col="green")
				lines(printorder[ord],mR[ord],col="red")
				lines(printorder[ord],mG[ord],col="green")
			}
			mR <- mR - mean(mR,na.rm=TRUE)
			mG <- mG - mean(mG,na.rm=TRUE)
		} else {
			plate.m <- tapply(Af,plate,hubermu)
			m <- Af
			for (i in 1:max(plate)) m[plate==i] <- plate.m[i]
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m,na.rm=TRUE)
		}
	} else {
		# Smooth correction for time order
		if(separate.channels) {
			mR <- fitted(loess(Rf~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			mG <- fitted(loess(Gf~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Rf,xlab="Print Order",ylab="Log Intensity",type="n")
				points(printorder,Rf,pch=".",col="red")
				points(printorder,Gf,pch=".",col="green")
				lines(printorder[ord],mR[ord],col="red")
				lines(printorder[ord],mG[ord],col="green")
			}
			mR <- mR - mean(mR,na.rm=TRUE)
			mG <- mG - mean(mG,na.rm=TRUE)
		} else {
			m <- fitted(loess(Af~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m,na.rm=TRUE)
		}
	}
	list(R=2^(Rf-mR),G=2^(Gf-mG),R.trend=mR,G.trend=mG)
}

plotPrintorder <- function(object,layout,start="topleft",slide=1,method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order.
#	Gordon Smyth
#	9 Apr 2002.  Last revised 18 June 2003.

	if(is.null(object$R) || is.null(object$G)) stop("List must contain components R and G")
	G <- object$G[,slide]
	R <- object$R[,slide]
	if(length(R) != length(G)) stop("R and G must have same length")
	start <- match.arg(start,c("topleft","topright"))
	po <- printorder(layout,start=start)$printorder
	invisible(normalizeForPrintorder.rg(R=R,G=G,printorder=po,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size,plot=TRUE))
}

#  BETWEEN ARRAY NORMALIZATION

setGeneric("normalizeBetweenArrays", function(object,method="scale") standardGeneric("normalizeBetweenArrays")) 

setMethod("normalizeBetweenArrays", "matrix", definition=
function(object, method="scale") {
#	Normalize between arrays - method for matrix
#	Gordon Smyth
#	12 Apri 2003.  Last revised 15 June 2003.

	choices <- c("none","scale","quantile")
	method <- choices[pmatch(method,choices)]
	if(is.na(method)) {
		warning("normalization method not recognized, defaulting to \"none\"")
		method <- "none"
	}
	switch(method,
		none = object,
		scale = normalizeMedianDeviations(object),
		quantile = normalizeQuantiles(object)
	)
})

setMethod("normalizeBetweenArrays", "list", definition=
function(object, method="scale") {
#	Normalize between arrays - method for list
#	Gordon Smyth
#	23 Apri 2003

#	Try to convert list to MAList
	normalizeBetweenArrays(new("MAList",unclass(object)), method=method)
})

setMethod("normalizeBetweenArrays", "MAList", definition=
function(object, method="scale") {
#	Normalize between arrays - method for MAList
#	Gordon Smyth
#	12 Apri 2003.  Last revised 15 June 2003.

	choices <- c("none","scale","quantile")
	method <- choices[pmatch(method,choices)]
	if(is.na(method)) warning("normalization method not recognized, defaulting to \"none\"")
	switch(method,
		scale = {
			object$M <- normalizeMedianDeviations(object$M)
			object$A <- normalizeMedians(object$A)
		},
		quantile = {
			narrays <- NCOL(object$M)
			Z <- normalizeQuantiles(cbind(object$A+object$M/2,object$A-object$M/2))
			R <- Z[,1:narrays]
			G <- Z[,narrays+(1:narrays)]
			object$M <- R-G
			object$A <- (R+G)/2
		})
	object
})

normalizeQuantiles <- function(A, ties=FALSE) {
#	Make all the columns of a matrix have the same quantiles, allowing for missing values.
#	Gordon Smyth
#	25 June 2002.  Last revised 5 June 2003.

	n <- dim(A)
	if(is.null(n)) return(A)
	if(n[2]==1) return(A)
	O <- S <- array(,n)
	if(ties) R <- O
	nobs <- rep(n[1],n[2])
	i <- (0:(n[1]-1))/(n[1]-1)
	for (j in 1:n[2]) {
		Si <- sort(A[,j], method="quick", index.return=TRUE)
		if(ties) R[,j] <- rank(A[,j])
		nobsj <- length(Si$x)
		if(nobsj < n[1]) {
			nobs[j] <- nobsj
			isna <- is.na(A[,j])
			S[,j] <- approx((0:(nobsj-1))/(nobsj-1), Si$x, i, ties="ordered")$y
			O[!isna,j] <- ((1:n[1])[!isna])[Si$ix]
		} else {
			S[,j] <- Si$x
			O[,j] <- Si$ix
		}
	}
	m <- rowMeans(S)
	for (j in 1:n[2]) {
		if(nobs[j] < n[1]) {
			isna <- is.na(A[,j])
			if(ties)
				A[!isna,j] <- approx(i, m, (R[!isna,j]-1)/(nobs[j]-1), ties="ordered")$y
			else
				A[O[!isna,j],j] <- approx(i, m, (0:(nobs[j]-1))/(nobs[j]-1), ties="ordered")$y
		} else {
			if(ties)
				A[,j] <- approx(i, m, (R[,j]-1)/(n[1]-1), ties="ordered")$y
			else
				A[O[,j],j] <- m
		}
	}
	A
}

normalizeMedianDeviations <- function(x) 
{
#	Normalize columns of a matrix to have the same median absolute value
#	Gordon Smyth
#	14 Mar 2002.  Last revised 12 Apr 2003.

	narrays <- NCOL(x)
	if(narrays==1) return(x)
	medabs <- function(x) median(abs(as.numeric(x[!is.na(x)])))
	xmat.mav <- apply(x, 2, medabs)
	denom <- (prod(xmat.mav))^(1/narrays)
	si <- xmat.mav/denom
	t(t(x)/si)
}

normalizeMedians <- function(x) 
{
#	Normalization columns of a matrix to have the same median value
#	Gordon Smyth
#	12 April 2003

	narrays <- NCOL(x)
	if(narrays==1) return(x)
	a.med <- apply(x, 2, median, na.rm=TRUE)
	a.med <- a.med / (prod(a.med))^(1/narrays)
	t(t(x)/a.med)
}
