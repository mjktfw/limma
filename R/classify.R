#  CLASSIFICATION

classifyTests <- function(tstat,cor.matrix=NULL,design=NULL,contrasts=diag(ncol(design)),df=Inf,p.value=0.01) {
#	Use F-tests to classify vectors of t-test statistics into outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 3 July 2003.

#	Method intended for MAList objects but accept unclassed lists as well
	if(is.list(tstat)) {
		if(is.null(tstat$t)) stop("tstat cannot be extracted from object")
		if(missing(design) && !is.null(tstat$design)) design <- tstat$design
		if(missing(contrasts) && !is.null(tstat$contrasts)) contrasts <- tstat$contrasts
		if(missing(df) && !is.null(tstat$df.prior) && !is.null(tstat$df.residual)) df <- tstat$df.prior+tstat$df.residual
		tstat <- tstat$t
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
	if(!is.null(design)) {
		design <- as.matrix(design)
		R <- chol(crossprod(design))
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

classifyTests43 <- function(tstat,t1=4,t2=3) {
#	Simple classification of vectors of t-test statistics
#	Gordon Smyth
#	1 July 2003.

	if(is.list(tstat)) tstat <- tstat$t
	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	apply(tstat,1,function(x) any(abs(x)>t1,na.rm=TRUE)) * sign(tstat)*(abs(tstat)>t2)
}
