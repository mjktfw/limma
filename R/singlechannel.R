#  SINGLE CHANNEL ANALYSIS

lmscFit <- function(object,design,correlation)
#	Fit single channel linear model for each gene to a series of microarrays
#	allowing for known correlation between the channels on each spot.
#	Gordon Smyth
#	14 March 2004.  Last modified 19 April 2004.
{
#	Check input
	M <- as.matrix(object$M)
	A <- as.matrix(object$A)
	if(is.null(M) || is.null(A)) stop("object must have components M and A")
	dimM <- dim(M)
	dimA <- dim(A)
	if(any(dimM != dimA)) stop("dimensions of M and A don't match")
	if(!all(is.finite(M)) || !all(is.finite(A))) stop("Missing or infinite values found in M or A")
	if(missing(design)) stop("design matrix must be specified")
	narrays <- dimM[2]
	ny <- 2*narrays
	design <- as.matrix(design)
	if(nrow(design) != ny) stop("The number of rows of the design matrix should match the number of channel intensities, i.e., twice the number of arrays")
	if(missing(correlation)) stop("intra-spot correlation must be specified")
	if(abs(correlation) >= 1) stop("correlation must be strictly between -1 and 1")

#	Dimensions
	nbeta <- ncol(design)
	coef.names <- colnames(design)
	ngenes <- dimM[1]

#	Main computation
	sdM <- sqrt(2*(1-correlation))
	sdA <- sqrt((1+correlation)/2)
	y <- rbind(t(M)/sdM, t(A)/sdA)
	designM <- (diag(narrays) %x% matrix(c(-1,1),1,2)) %*% design
	designA <- (diag(narrays) %x% matrix(c(0.5,0.5),1,2)) %*% design
	X <- rbind(designM/sdM, designA/sdA)

#	In general it may be necessary to allow for quality weights, this call does not
	fit <- lm.fit(X,y)
	fit$sigma <- sqrt(colSums(fit$effects[(fit$rank+1):ny,]^2) / fit$df.residual)
	fit$fitted.values <- fit$residuals <- fit$effects <- NULL
#	if(variance.smooth) fit$s2 <- squeezeVar(fit$s2, fit$df.residual)
	fit$coefficients <- t(fit$coefficients)
	stdev.unscaled <- sqrt(diag(chol2inv(fit$qr$qr)))
	fit$stdev.unscaled <- matrix(stdev.unscaled,ngenes,nbeta,byrow=TRUE)
	fit$df.residual <- rep.int(fit$df.residual,ngenes)
	dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
	fit$design <- design
	fit$correlation <- correlation
	fit$genes <- object$genes
	fit$Amean <- rowMeans(A,na.rm=TRUE)
	new("MArrayLM",fit)
}

intraspotCorrelation <- function(object,design,trim=0.15)
#	Estimate intra-spot correlation between channels for two channel data
#	Gordon Smyth
#	19 April 2004.
{
#	Check input
	M <- as.matrix(object$M)
	A <- as.matrix(object$A)
	if(is.null(M) || is.null(A)) stop("object should have components M and A")
	dimM <- dim(M)
	dimA <- dim(A)
	if(any(dimM != dimA)) stop("dimensions of M and A don't match")
	if(!all(is.finite(M)) || !all(is.finite(A))) stop("Missing or infinite values found in M or A")
	if(missing(design)) stop("design matrix must be specified")
	ngenes <- dimM[1]
	narrays <- dimM[2]
	ny <- 2*narrays
	design <- as.matrix(design)
	if(nrow(design) != ny) stop("The number of rows of the design matrix should match the number of channel intensities, i.e., twice the number of arrays")

#	Fit heteroscedastic regression for each gene
	Ident <- diag(narrays)
	designM <- (Ident %x% matrix(c(-1,1),1,2)) %*% design
	designA <- (Ident %x% matrix(c(0.5,0.5),1,2)) %*% design
	X <- rbind(designM, designA)
	Z <- diag(2) %x% rep(1,narrays)
	if(!require(statmod)) stop("statmod package required but is not available")
	arho <- rep(0,ngenes)
	degfre <- matrix(0,ngenes,2,dimnames=list(rownames(M),c("df.M","df.A")))
	for (i in 1:ngenes) {
		y <- c(M[i,],A[i,])
		fit <- remlscore(y,X,Z)
		arho[i] <- 0.5*(fit$gamma[2]-fit$gamma[1])
		degfre[i,] <- crossprod(Z,1-fit$h)
	}
	arho <- arho+log(2)
	arhobias <- digamma(degfre[,1]/2)-log(degfre[,1]/2)-digamma(degfre[,2]/2)+log(degfre[,2]/2)
	list(consensus.correlation=tanh(mean(arho-arhobias,trim=trim)), all.correlations=arho, df=degfre)
}

array2channel <- function(targets,channels=c(1,2),channelwise.columns=list(Target=c("Cy3","Cy5")))
#	Convert data.frame with one row for each two-color array,
#	into data.frame with one row for each channel
#	Gordon Smyth
#	16 March 2004.  Last modified 13 May 2004.
{
	targets <- as.data.frame(targets)
	if(!min(dim(targets))) return(targets)
	lcc <- length(channelwise.columns)
	if(lcc) {
		out <- channelwise.columns
		cheaders <- names(channelwise.columns)
		for (i in 1:lcc) {
			aheaders <- channelwise.columns[[i]]
			if(all(aheaders %in% names(targets))) {
				out[[i]] <- as.vector(as.matrix((targets[,aheaders])))
				targets[[ aheaders[1] ]] <- NULL
				targets[[ aheaders[2] ]] <- NULL
			} else {
				out[[ cheaders[i] ]] <- NULL
			}
		}
	}
	narrays <- nrow(targets)
	channel.col <- rep(channels,c(narrays,narrays))
	out <- cbind(Channel=channel.col,rbind(targets,targets),out)
	row.names(out) <- paste(row.names(targets),channel.col,sep=".")
	o <- as.vector(t(matrix(1:(2*narrays),narrays,2)))
	out[o,]
}

