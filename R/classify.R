#  CLASSIFICATION

classifyTests <- function(tstat,cor.matrix=diag(ntests),df.residual=Inf,p.value=0.01) {
#	Classify a series of vectors of t-test statistics into vector outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 11 June 2003.

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
	if(length(qF)==1) qF <- rep(qF,ngenes) 
	result <- matrix(0,ngenes,ntests,dimnames=dimnames(tstat))
	for (i in 1:ngenes) {
		x <- tstat[i,]
		if( crossprod(crossprod(Q,x)) > qF[i] ) {
			ord <- order(abs(x),decreasing=TRUE)
			result[i,ord[1]] <- sign(x[ord[1]])
			for (j in 2:ntests) {
				bigger <- ord[1:(j-1)]
				x[bigger] <- sign(x[bigger]) * abs(x[ord[j]])
				if( crossprod(crossprod(Q,x)) > qF[i] )
					result[i,ord[j]] <- sign(x[ord[j]])
				else
					break
			}
		}
	}
	result
}
