#  CLASSIFICATION

classifyTests <- function(tstat,cor.matrix=diag(ntests),df.residual=Inf,p.value=0.01) {
#	Classify a series of vectors of t-test statistics into vector outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 9 June 2003.

	if(any(is.na(tstat))) stop("missing values not allowed")
	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	ngenes <- nrow(tstat)
	ntests <- ncol(tstat)
	if(ntests == 1) {
		qT <- qt(p.value/2,df.residual,lower.tail=FALSE)
		return( sign(tstat) * (abs(tstat) > qT) )
	}

#	cor.matrix is estimated correlation matrix of the coefficients
#	and also the estimated covariance matrix of the t-statistics
	E <- eigen(cor.matrix,symmetric=TRUE)
	r <- sum(E$values/E$values[1] > 1e-8)
	Q <- E$vectors[,1:r]

	qF <- r * qf(p.value, r, df.residual, lower.tail=FALSE)
	if(length(pF)==1) pF <- rep(pF,ngenes) 
	result <- tstat
	for (i in 1:ngenes) {
		x <- tstat[i,]
		if( crossprod(crossprod(Q,x)) > qF ) {
			ord <- order(abs(x),decreasing=TRUE)
			result[i,ord[1]] <- sign(x[ord[1]])
			for (i in 2:p) {
				bigger <- ord[1:(i-1)]
				x[bigger] <- sign(x[bigger]) * abs(x[ord[i]])
				if( crossprod(crossprod(Q,x)) > qF )
					result[ord[i]] <- sign(x[ord[i]])
				else
					break
			}
		}
	}
	result
}

classifyTests.matrix <- function(tstat,cor.matrix=diag(ntests),df.residual=Inf,p.value=0.01) {
#	Classify a series of vectors of t-test statistics into vector outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 9 June 2003.

	if(any(is.na(tstat))) stop("missing values not allowed")
	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	ngenes <- nrow(tstat)
	ntests <- ncol(tstat)

#	cor.matrix is estimated correlation matrix of the coefficients
#	and also the estimated covariance matrix of the t-statistics
	E <- eigen(cor.matrix,symmetric=TRUE)
	r <- sum(E$values/E$values[1] > 1e-8)
	Q <- E$vectors[,1:r]

	pF <- r * qf(p.value, r, df.residual, lower.tail=FALSE)
	if(length(pF)==1) pF <- rep(pF,ngenes) 
	result <- tstat
	for (i in 1:ngenes)
		result[i,] <- classifyTests.vector(tstat[i,], Q, pF[i])
	result
}

classifyTests.vector <- function(tstat, Q=diag(p), pF=ncol(Q)*4) {
#	Classify a vector of t-test statistics into a vector outcome
#	Gordon Smyth
#	20 Mar 2003.  Last modified 9 June 2003.

	p <- length(tstat)
	result <- rep(0,p)
	if( crossprod(crossprod(Q,tstat)) > pF ) {
		ord <- order(abs(tstat),decreasing=TRUE)
		result[ord[1]] <- sign(tstat[ord[1]])
		if(p > 1) for (i in 2:p) {
			x <- tstat
			bigger <- ord[1:(i-1)]
			x[bigger] <- sign(x[bigger]) * abs(x[ord[i]])
			if( crossprod(crossprod(Q,x)) > pF )
				result[ord[i]] <- sign(tstat[ord[i]])
			else
				break
		}
	}
	result
}
