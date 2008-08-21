##  GENESET.R

geneSetTest <- function(selected,statistics,alternative="mixed",type="auto",ranks.only=TRUE,nsim=10000)
#	Gene set test using either Wilcox test or simulation.
#	Gordon Smyth
#	3 September 2004. Last modified 21 July 2006.
{
	alternative <- match.arg(alternative,c("mixed","either","down","up","less","greater"))
	if(alternative=="two.sided") alternative <- "either"
	if(alternative=="less") alternative <- "down"
	if(alternative=="greater") alternative <- "up"
	type <- match.arg(tolower(type),c("auto","t","f"))
	allsamesign <- all(statistics >= 0) || all(statistics <= 0)
	if(type=="auto") {
		if(allsamesign)
			type <- "f"
		else
			type <- "t"
	}
	if(type=="f" & alternative!="mixed") stop("Only alternative=\"mixed\" is possible with F-like statistics.")
	if(alternative=="mixed") statistics <- abs(statistics)
	if(alternative=="down") {
		statistics <- -statistics
		alternative <- "up"
	}
	if(ranks.only) {
#		The test statistic is the mean rank of the selected statistics
#		and the p-value is obtained explicitly from the Wilcox test
		if(alternative=="either")
			wilc.alt <- "two.sided"
		else
			wilc.alt <- "greater"
		x <- y <- NULL
		if(is.logical(selected)) {
			x <- statistics[selected]
			y <- statistics[!selected]
		}
		if(is.numeric(selected)) {
			x <- statistics[selected]
			y <- statistics[-selected]
		}
		if(is.character(selected)) {
			nam <- names(statistics)
			if(is.null(nam)) stop("selected is character but elements of statistics are not named")
			selected <- is.element(nam,selected)
			x <- statistics[selected]
			y <- statistics[!selected]
		}
		return(wilcox.test(x,y,alternative=wilc.alt,conf.int=FALSE)$p.value)
	} else {
#		The test statistic is the mean of the selected statistics
#		and the p-value is obtained by random permutation
		ssel <- statistics[selected]
		ssel <- ssel[!is.na(ssel)]
		nsel <- length(ssel)
		if(nsel==0) return(1)
		stat <- statistics[!is.na(statistics)]
		msel <- mean(ssel)
		if(alternative=="either")
			posstat <- abs
		else
			posstat <- function(x) x
		msel <- posstat(msel)
		ntail <- 0
		for (i in 1:nsim) if(posstat(mean(sample(stat,nsel))) >= msel) ntail <- ntail+1
		return(ntail/nsim)
	}
}


roast <- function(iset=NULL,y,design,contrast=ncol(design),gene.weights=NULL,array.weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,nrot=1000)
# rotation gene set testing for linear models
# Gordon Smyth and Di Wu
# 24 Apr 2008. Last revised 7 Aug 2008.
{
	if(is.null(iset)) iset <- rep(TRUE,nrow(y))
	y <- as.matrix(y)
	design <- as.matrix(design)
	if(!is.null(df.prior) && df.prior<0) stop("df.prior must be non-negative")
	
	p <- ncol(design)
	p0 <- p-1
	n <- ncol(y)
	d <- n-p

	if(length(contrast) == 1) {
		u <- rep.int(0,p)
		u[contrast] <- 1
		contrast <- u
	}
	if(length(contrast) != p) stop("length of contrast must match column dimension of design")
	if(all(contrast==0)) stop("contrast all zero")
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")

	if(!is.null(array.weights)) {
		if(any(array.weights <= 0)) stop("array.weights must be positive")
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		design <- design*sqrt(array.weights)
		y <- t(t(y)*sqrt(array.weights))
	}

	if(!is.null(block)) {
		block <- as.vector(block)
		if (length(block) != n) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,n,nblocks) == matrix(ub,n,nblocks,byrow = TRUE)
		cormatrix <- Z %*% (correlation * t(Z))
		diag(cormatrix) <- 1
		R <- chol(cormatrix)
		y <- t(backsolve(R, t(y), transpose = TRUE))
		design <- backsolve(R, design, transpose = TRUE)
 	}

	qr <- qr(contrast)
	Q <- qr.Q(qr,complete=TRUE)
	sign1 <- sign(qr$qr[1,1])
	Q <- cbind(Q[,-1],Q[,1])
	X <- design %*% Q
	qr <- qr(X)
	sign2 <- sign(qr$qr[p,p])
	signc <- sign1 * sign2

	if(is.null(var.prior) || is.null(df.prior)) {
#		Fit model to all genes
		effects <- qr.qty(qr,t(y))
#		Estimate global parameters s0 and d0
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		sv <- squeezeVar(s2,df=d)
		d0 <- sv$df.prior
		s02 <- sv$var.prior
		effects <- effects[,iset,drop=FALSE]
		sd.post <- sqrt(sv$var.post[iset])
	} else {
		y <- y[iset,,drop=FALSE]
		effects <- qr.qty(qr,t(y))
		d0 <- df.prior
		s02 <- var.prior
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(is.finite(d0))
			sd.post <- sqrt( (d0*s02+d*s2)/(d0+d) )
		else
			sd.post <- sqrt(s02)
	}

#	From here, all results are for set only
	nset <- ncol(effects)
	Y <- effects[-(1:p0),,drop=FALSE]
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post

	statobs <- p <- rep(0,4)
	names(statobs) <- names(p) <- c("mixed","up","down","either")
	stati <- array(0,c(nrot,4),dimnames=list(NULL,names(p)))

#	Convert to z-scores
	modt <- zscoreT(modt,df=d0+d)
	if(!is.null(gene.weights)) {
		if(length(gene.weights) != nset) stop("length of gene.weights disagrees with size of set")
		w <- sign(gene.weights)*sqrt(abs(gene.weights))
		modt <- w*modt
	}

#	Observed statistics
	modt2 <- modt^2
	mp <- sum(modt2[modt>0])/nset
	mn <- sum(modt2[modt<0])/nset
	statobs["mixed"] <- mean(modt2) 
	statobs["up"] <- mp
	statobs["down"] <- mn
	statobs["either"] <- max(c(mp,mn))

#	Random rotations
	R <- matrix(rnorm(nrot*(d+1)),nrot,d+1)
	R <- R/sqrt(rowSums(R^2))
	Br <- R %*% Y
	s2r <- (matrix(YY,nrot,nset,byrow=TRUE)-Br^2)/d
	if(is.finite(d0))
		sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
	else
		sdr.post <- sqrt(s02)
	modtr <- signc*Br/sdr.post
	modtr <- zscoreT(modtr,df=d0+d)
	if(!is.null(gene.weights)) modtr <- t(w*t(modtr))

	mp <- rowMeans(pmax(modtr,0)^2)
	mn <- rowMeans(pmax(-modtr,0)^2)
	stati[,"mixed"] <- rowMeans(modtr^2)
	stati[,"up"] <- mp
	stati[,"down"] <- mn
	stati[,"either"] <- pmax(mp,mn)

#	p-values
	p <- rowMeans(t(stati) >= statobs)

#	Output
	r1 <- mean(modt > 1)
	r2 <- mean(modt < -1)
	out <- data.frame(Z=sqrt(statobs),Active=c(r1+r2,r1,r2,max(r1,r2)),P.Value=p)
	row.names(out) <- c("mixed","up","down","either")
	out
}
