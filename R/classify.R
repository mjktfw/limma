#  CLASSIFICATION

classifyTests <- function(tstat,cor.matrix=NULL,design=NULL,contrasts=NULL,df=Inf,p.value=0.01) {
#	Classify a series of vectors of t-test statistics into vector outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 19 June 2003.

	if(is(tstat,"MArrayLM")) {
		cor.matrix <- NULL
		design <- tstat@design
		contrasts <- tstat@contrasts
		df <- tstat@df.prior+tstat@df.residual
		tstat <- tstat@tstat
	}

	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	ngenes <- nrow(tstat)
	ntests <- ncol(tstat)
	if(ntests == 1) {
		qT <- qt(p.value/2, df, lower.tail=FALSE)
		return( sign(tstat) * (abs(tstat) > qT) )
	}

#	cor.matrix is estimated correlation matrix of the coefficients
#	and also the estimated covariance matrix of the t-statistics
	if(!is.null(cor.matrix) && !is.null(design)) stop("Cannot specify both cor.matrix and design")
	if(is.null(cor.matrix) && !is.null(design)) {
		R <- chol(crossprod(design))
		if(length(contrasts)==0) contrasts <- diag(ncol(design))
		cor.matrix <- crossprod(backsolve(R,contrasts,transpose=TRUE))
		d <- sqrt(diag(cor.matrix))
		cor.matrix <- cor.matrix / (d %*% t(d))
	}
	if(is.null(cor.matrix)) {
		r <- ntests
		Q <- diag(r)
	} else {
		E <- eigen(cor.matrix,symmetric=TRUE)
		r <- sum(E$values/E$values[1] > 1e-8)
		Q <- matvec( E$vectors[,1:r], 1/sqrt(E$values[1:r]))
	}

	qF <- r * qf(p.value, r, df, lower.tail=FALSE)
	if(length(qF)==1) qF <- rep(qF,ngenes) 
	result <- matrix(0,ngenes,ntests,dimnames=dimnames(tstat))
	if(is.null(colnames(tstat)) && !is.null(colnames(contrasts))) colnames(result) <- colnames(contrasts)
	for (i in 1:ngenes) {
		x <- tstat[i,]
		if(any(is.na(x)))
			result[i,] <- NA
		else
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
