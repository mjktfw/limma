#  WITHIN ARRAY NORMALIZATION

require(modreg)

MA.RG <- function(object) {
#	Convert RGList to MAList
#	Gordon Smyth
#	2 March 2003.  Last revised 26 April 2003.

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
	object$M <- R-G
	object$A <- (R+G)/2
	new("MAList",unclass(object))
}

if(!isGeneric("normalizeWithinArrays"))
	setGeneric("normalizeWithinArrays", function(object,...) standardGeneric("normalizeWithinArrays")) 

setMethod("normalizeWithinArrays", "list", definition=
function(object, ...) {
	normalizeWithinArrays(new("RGList",object), ...)
})

setMethod("normalizeWithinArrays", "RGList", definition=
function(object,layout,method="printtiploess",weights=object$weights,span=0.3,iterations=4,controlspots=NULL) {
#	Sub-array loess normalization
#	Gordon Smyth
#	2 March 2003.  Last revised 25 April 2003.

	MA <- MA.RG(object)
	choices <- c("median","loess","printtiploess","composite","plate","printorder")
	method <- choices[pmatch(method,choices)]
	if(is.na(method)) return(MA)
	ngenes <- NROW(MA$M)
	narrays <- NCOL(MA$M)
	switch(method,
		median = {
			for (j in 1:narrays) MA$M[,j] <- MA$M[,j] - median(MA$M[,j],na.rm=TRUE)
		},
		loess = {
			for (j in 1:narrays) {
				y <- MA$M[,j]
				x <- MA$A[,j]
				w <- weights[,j]
				MA$M[,j] <- residuals(loess(y~x,weights=w,span=0.75*span,na.action=na.exclude,degree=1,family="symmetric",trace.hat="approximate",iterations=iterations,surface="direct"))
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
					y <- MA$M[spots,j]
					x <- MA$A[spots,j]
					w <- weights[spots,j]
					MA$M[spots,j] <- residuals(loess(y~x,weights=w,span=span,na.action=na.exclude,degree=1,family="symmetric",trace.hat="approximate",iterations=iterations,surface="direct"))
					spots <- spots + nspots
				}
			}
		},
		composite = {
			for (j in 1:narrays) {
				y <- MA$M[,j]
				x <- MA$A[,j]
				w <- weights[,j]
				fit <- loess(y~x,weights=w,span=span,subset=controlspots,na.action=na.exclude,degree=0,surface="direct",family="symmetric",trace.hat="approximate",iterations=iterations)
				global <- predict(fit,newdata=x)
				alpha <- (rank(x)-1) / sum(!is.na(x))
				ntips <- layout$ngrid.r * layout$ngrid.c
				nspots <- layout$nspot.r * layout$nspot.c
				spots <- 1:nspots
				for (tip in 1:ntips) {
					y <- MA$M[spots,j]
					x <- MA$A[spots,j]
					w <- weights[spots,j]
					local <- fitted(loess(y~x,weights=w,span=span,na.action=na.exclude,degree=1,family="symmetric",trace.hat="approximate",iterations=iterations,surface="direct"))
					MA$M[spots,j] <- MA$M[spots,j] - alpha[spots]*global[spots]-(1-alpha[spots])*local
					spots <- spots + nspots
				}
			}
		}
	)
	MA
})

setGeneric("normalizeForPrintorder",function(object,...) standardGeneric("normalizeForPrintorder"))

setMethod("normalizeForPrintorder", "RGList",
function(object, ...) {
#	Pre-normalize the foreground intensities for print order - method for RGList
#	Gordon Smyth
#	11 Mar 2002.  Last revised 24 April 2003.

	nslides <- NCOL(object$R)
	for (i in 1:nslides) {
		RG <- normalizeForPrintorder(object$R[,i],object$G[,i], ...)
		object$R[,i] <- RG$R
		object$G[,i] <- RG$G
	}
	object
})

setMethod("normalizeForPrintorder", "ANY",
function(object, ...) normalizeForPrintorder(R=object, ...)
)

setMethod("normalizeForPrintorder", "missing",
function(R,G,layout,method="loess",separate.channels=FALSE,span=0.1,plate.size=32,plot=FALSE) {
#	Pre-normalize the foreground intensities for print order, given R and G for a single array.
#	Gordon Smyth
#	8 Mar 2002.  Last revised 24 Apr 2003.

	ngrid.r <- layout$ngrid.r
	ngrid.c <- layout$ngrid.c
	nspot.r <- layout$nspot.r
	nspot.c <- layout$nspot.c
	spot.c <- rep(1:nspot.c,ngrid.r*ngrid.c*nspot.r)
	spot.r <- rep(rep(1:nspot.r,rep(nspot.c,nspot.r)),ngrid.r*ngrid.c)
	printorder <- nspot.c*(spot.r-1)+nspot.c+1-spot.c
	if(plot) ord <- order(printorder)
	Rf <- log(R,2)
	Gf <- log(G,2)
	if(!separate.channels) Af <- (Rf+Gf)/2
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
			mR <- mR - mean(mR)
			mG <- mG - mean(mG)
		} else {
			plate.m <- tapply(Af,plate,hubermu)
			m <- Af
			for (i in 1:max(plate)) m[plate==i] <- plate.m[i]
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m)
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
			mR <- mR - mean(mR)
			mG <- mG - mean(mG)
		} else {
			m <- fitted(loess(Af~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m)
		}
	}
	list(R=2^(Rf-mR),G=2^(Gf-mG),printorder=printorder,R.trend=mR,G.trend=mG)
})

plotPrintorder <- function(R,G=NULL,layout,slide=1,method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order.
#	Input can be an RG list for a series of arrays or an R,G pair for a single array.
#	Gordon Smyth
#	9 Apr 2002.

	if(is.list(R)) {
		if(is.null(R$R) || is.null(R$G)) stop("List must contain components R and G")
		G <- R$G[,slide]
		R <- R$R[,slide]
	}
	if(length(R) != length(G)) stop("R and G must have same length")
	invisible(normalizeForPrintorder(R=R,G=G,layout=layout,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size,plot=TRUE))
}

#  BETWEEN ARRAY NORMALIZATION

if(!isGeneric("normalizeBetweenArrays"))
	setGeneric("normalizeBetweenArrays", function(object,...) standardGeneric("normalizeBetweenArrays")) 

setMethod("normalizeBetweenArrays", "matrix", definition=
function(object, method="scale") {
#	Normalize between arrays - method for matrix
#	Gordon Smyth
#	12 Apri 2003

	choices <- c("none","scale","quantile")
	method <- choices[pmatch(method,choices)]
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
	normalizeBetweenArrays(new("MAList",object), method=method)
})

setMethod("normalizeBetweenArrays", "MAList", definition=
function(object, method="scale") {
#	Normalize between arrays - method for MA-list
#	Gordon Smyth
#	12 Apri 2003

	choices <- c("none","scale","quantile")
	method <- choices[pmatch(method,choices)]
	switch(method,
		scale = {
			object$M <- normalizeMedianDeviations(object$M)
			object$A <- normalizeMedians(object$A)},
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

normalizeQuantiles <- function(A) {
#	Make all the columns of a matrix have the same quantiles, allowing for missing values.
#	Gordon Smyth
#	25 June 2002.  Last revised 12 Mar 2003.

	n <- dim(A)
	if(is.null(n)) return(A)
	if(n[2]==1) return(A)
	O <- S <- rep(NA,n[1]*n[2])
	dim(O) <- dim(S) <- n
	for (j in 1:n[2]) {
		Si <- sort(A[,j],index=TRUE)
		isna <- is.na(A[,j])
		if(any(isna)) {
			nobs <- length(Si$x)
			S[,j] <- approx((0:(nobs-1))/(nobs-1), Si$x, (0:(n[1]-1))/(n[1]-1))$y
			O[!isna,j] <- ((1:n[1])[!isna])[Si$ix]
		} else {
			S[,j] <- Si$x
			O[,j] <- Si$ix
		}
	}
	m <- rowMeans(S)
	for (j in 1:n[2]) {
		isna <- is.na(A[,j])
		if(any(isna)) {
			nobs <- sum(!isna)
			A[O[!isna,j],j] <- approx((0:(n[1]-1))/(n[1]-1), m, (0:(nobs-1))/(nobs-1))$y
		} else
			A[O[,j],j] <- m 
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
