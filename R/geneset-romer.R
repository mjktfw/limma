##  ROMER.R

symbols2indices <- function(gmtl.official, symbol)
# Make a list of gene sets symbols into a list of gene sets indices used to create input for romer
# Gordon Smyth and Yifang Hu
# 25 March 2009.
{
	iset<-list()
	j<-1
	for(i in 1:length(gmtl.official))
	{
		idx<-match(gmtl.official[[i]],symbol)
		if(any(is.na(idx)=="FALSE"))
		{
			iset[[j]] <-idx[!is.na(idx)]
			names(iset)[j]<-names(gmtl.official)[i]
			j<-j+1
		}
	}

	iset
}

romer2 <- function(iset,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,nrot=9999)
# rotation-mean50-rank version of GSEA (gene set enrichment analysis) for linear models
# Gordon Smyth and Yifang Hu
# 27 March 2009.  Last modified 15 Sep 2009.
{
#	Check input arguments
	if(!is.list(iset)) iset <- list(set=iset)
	nset<-length(iset)
	y <- as.matrix(y)
	design <- as.matrix(design)
	ngenes<-nrow(y)	
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

#	Divide out array weights, if they exist
	if(!is.null(array.weights)) {
		if(any(array.weights <= 0)) stop("array.weights must be positive")
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		design <- design*sqrt(array.weights)
		y <- t(t(y)*sqrt(array.weights))
	}

#	Divide out correlation structure, it is exists
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

#	Reform design matrix so that contrast of interest is last column
	qr <- qr(contrast)
	Q <- qr.Q(qr,complete=TRUE)
	sign1 <- sign(qr$qr[1,1])
	Q <- cbind(Q[,-1],Q[,1])
	X <- design %*% Q
	qr <- qr(X)
	sign2 <- sign(qr$qr[p,p])
	signc <- sign1 * sign2

# 	Fit model to all genes
	effects <- qr.qty(qr,t(y))

# 	Estimate global hyper-parameters s0 and d0
	s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
	sv <- squeezeVar(s2,df=d)
	d0 <- sv$df.prior
	s02 <- sv$var.prior
	sd.post <- sqrt(sv$var.post)

#	t-statistics and effects
	Y <- effects[-(1:p0),,drop=FALSE]
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post
	
#	Observed rankings for each set
	s.r <- rank(modt)

	s.rank.mixed<-rep(0,nset)
	s.rank.up<-rep(0,nset)
	s.rank.down<-rep(0,nset)
	s.rank.either<-rep(0,nset)

	modt.abs<-abs(modt)
	s.abs.r <-rank(modt.abs)

	m<-unlist(lapply(iset,length))
	m<-floor((m+1)/2)

	for(i in 1:nset)
	{	
		mh<-.meanHalf(s.r[iset[[i]]],m[i])
		s.rank.up[i] <-mh[2]	
		s.rank.down[i]<-mh[1]
  		s.rank.either[i]<- max(abs(mh-(ngenes+1)/2))
  		s.rank.mixed[i]<-.meanHalf(s.abs.r[iset[[i]]],m[i])[2]
	}	

#	Estimate hyper-parameters p0
	p.obs <- 2*pt(abs(modt),df=d0+d,lower.tail=FALSE)
	p0.obs <- 1-convest(p.obs) # proportion of DE probes
	
#	Estimate hyper-paremeter v0
	covmat <- chol2inv(qr$qr, size = qr$rank)
 	stdev.unscaled <- rep(sqrt(covmat[qr$rank,qr$rank]),ngenes)
	
	proportion<-p0.obs
	stdev.coef.lim <- c(0.1, 4)
	
	# get v0
	df.total <- rep(d,ngenes) + sv$df.prior
	var.prior.lim <- stdev.coef.lim^2/sv$var.prior
	var.prior <- tmixture.vector(modt, stdev.unscaled, df.total,proportion, var.prior.lim)
	if (any(is.na(var.prior))) {
		var.prior[is.na(var.prior)] <- 1/sv$var.prior
		warning("Estimation of var.prior failed - set to default value")
	}
	
#	Estimate posterior probability of DE
	r <- rep(1, ngenes) %o% var.prior
	r <- (stdev.unscaled^2 + r)/stdev.unscaled^2

	if (sv$df.prior > 10^6)
		kernel <- modt^2 * (1 - 1/r)/2
	else
		kernel <- (1 + df.total)/2 * log((modt^2 + df.total)/(modt^2/r + df.total))

	lods <- log(proportion/(1 - proportion)) - log(r)/2 + kernel
	pg <- exp(lods)/(1+exp(lods))

#	Shrink contrast to be like a residual
	Y[1,] <- Y[1,]*(sv$var.post/(sv$var.post+var.prior*pg))^(1/2)

#	Random rotations
	p.value <- matrix(rep(0,nset*4),nrow=nset,ncol=4)
	for(i in 1:nrot)
	{
		R <- matrix(rnorm((d+1)),1,d+1)
		R <- R/sqrt(rowSums(R^2))
		Br <- R %*% Y
		s2r <- (YY-Br^2)/d

		if(is.finite(d0))
			sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
		else
			sdr.post <- sqrt(s02)

		modtr <- signc*Br/sdr.post
		modtr.abs<-abs(modtr)
	
		s.r.2<-rank(modtr)
		s.abs.r.2<-rank(modtr.abs)
	
		for(j in 1:nset)
		{
			mh.2<-.meanHalf(s.r.2[iset[[j]]],m[j])
			
			s.rank.up.2 <-mh.2[2]	
			s.rank.down.2 <-mh.2[1]
  			s.rank.mixed.2 <-.meanHalf(s.abs.r.2[iset[[j]]],m[j])[2]
			s.rank.either.2 <- max(abs(mh.2-(ngenes+1)/2))
		
			if(s.rank.mixed.2>=s.rank.mixed[j]) p.value[j,1]<-p.value[j,1]+1
			if(s.rank.up.2>=s.rank.up[j]) p.value[j,2]<-p.value[j,2]+1
			if(s.rank.down.2<=s.rank.down[j]) p.value[j,3]<-p.value[j,3]+1
			if(s.rank.either.2>=s.rank.either[j]) p.value[j,4]<-p.value[j,4]+1
		}
	}	

	p.value <- (p.value+1)/(nrot+1)
	colnames(p.value)<-c("mixed","up","down","either")
	SetNames <- names(iset)
	if(is.null(SetNames))
		rownames(p.value) <- 1:nset
	else
		rownames(p.value) <- SetNames
	len.iset<-as.numeric(lapply(iset,length))
	cbind(NGenes=len.iset,p.value)
}

## Return means of top half and bottom half of the ranks for romer2
.meanHalf<- function(x,n)
#	Yifang Hu
{
	l<-length(x)
	a<-sort(x,partial=n)
	top<-mean(a[1:n])
	if((l%%2)==0) bottom<-mean(a[(n+1):l])
	if((l%%2)!=0) bottom<-mean(a[n:l])
	c(top,bottom)
}

romer <- function(iset,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,floor=FALSE,nrot=9999)
# rotation mean-rank version of GSEA (gene set enrichment analysis) for linear models
# Gordon Smyth and Yifang Hu
# 27 March 2009.  Last modified 5 Oct 2009.
{
#	Check input arguments
	if(!is.list(iset)) iset <- list(set=iset)
	nset<-length(iset)
	y <- as.matrix(y)
	design <- as.matrix(design)
	ngenes<-nrow(y)	
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

#	Divide out array weights, it they exist
	if(!is.null(array.weights)) {
		if(any(array.weights <= 0)) stop("array.weights must be positive")
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		design <- design*sqrt(array.weights)
		y <- t(t(y)*sqrt(array.weights))
	}

#	Divide out correlation structure, it is exists
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

#	Reform design matrix so that contrast of interest is last column
	qr <- qr(contrast)
	Q <- qr.Q(qr,complete=TRUE)
	sign1 <- sign(qr$qr[1,1])
	Q <- cbind(Q[,-1],Q[,1])
	X <- design %*% Q
	qr <- qr(X)
	sign2 <- sign(qr$qr[p,p])
	signc <- sign1 * sign2

# 	Fit model to all genes
	effects <- qr.qty(qr,t(y))

#	Sample variances
	s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)

# 	Estimate global hyper-parameters s0 and d0
	sv <- squeezeVar(s2,df=d)
	d0 <- sv$df.prior
	s02 <- sv$var.prior
	sd.post <- sqrt(sv$var.post)

#	t-statistics and effects (orthogonal residuals)
	Y <- effects[-(1:p0),,drop=FALSE]
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post
	
#	Estimate hyper-parameter p0
	p.obs <- 2*pt(abs(modt),df=d0+d,lower.tail=FALSE)
	p0.obs <- 1-convest(p.obs) # proportion of DE probes
	
#	Estimate hyper-paremeter v0
	covmat <- chol2inv(qr$qr, size = qr$rank)
 	stdev.unscaled <- rep(sqrt(covmat[qr$rank,qr$rank]),ngenes)
	proportion<-p0.obs
	stdev.coef.lim <- c(0.1, 4)
	df.total <- rep(d,ngenes) + sv$df.prior
	var.prior.lim <- stdev.coef.lim^2/sv$var.prior
	var.prior <- tmixture.vector(modt, stdev.unscaled, df.total,proportion, var.prior.lim)
	if (any(is.na(var.prior))) {
		var.prior[is.na(var.prior)] <- 1/sv$var.prior
		warning("Estimation of var.prior failed - set to default value")
	}
	
#	Estimate posterior probability of DE for each probe
	r <- rep(1, ngenes) %o% var.prior
	r <- (stdev.unscaled^2 + r)/stdev.unscaled^2
	if (sv$df.prior > 10^6)
		kernel <- modt^2 * (1 - 1/r)/2
	else
		kernel <- (1 + df.total)/2 * log((modt^2 + df.total)/(modt^2/r + df.total))
	lods <- log(proportion/(1 - proportion)) - log(r)/2 + kernel
	pg <- exp(lods)/(1+exp(lods))
	
#	Observed rankings for each set
	obs.ranks <- matrix(0,ngenes,3)
	if(floor) {
		obs.ranks[,1] <- rank(pmax(modt,0))
		obs.ranks[,2] <- rank(pmax(-modt,0))
		obs.ranks[,3] <- rank(pmax(abs(modt),1))
	} else {
		obs.ranks[,1] <- rank(modt)
		obs.ranks[,2] <- ngenes-obs.ranks[,1]+1
		obs.ranks[,3] <- rank(abs(modt))
	}
	obs.set.ranks <- matrix(0,nset,3)
	for(i in 1:nset) obs.set.ranks[i,] <- colMeans(obs.ranks[iset[[i]],,drop=FALSE])

#	Shrink contrast to be like a residual
	Y[1,] <- Y[1,]*(sv$var.post/(sv$var.post+var.prior*pg))^(1/2)

#	Random rotations to simulate null hypothesis
	rot.ranks <- obs.ranks
	p.value <- matrix(0,nrow=nset,ncol=3)
	for(i in 1:nrot)
	{
		R <- matrix(rnorm((d+1)),1,d+1)
		R <- R/sqrt(rowSums(R^2))
		Br <- R %*% Y
		s2r <- (YY-Br^2)/d

		if(is.finite(d0))
			sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
		else
			sdr.post <- sqrt(s02)

		modtr <- signc*Br/sdr.post
	
		if(floor) {
			rot.ranks[,1] <- rank(pmax(modtr,0))
			rot.ranks[,2] <- rank(pmax(-modtr,0))
			rot.ranks[,3] <- rank(pmax(abs(modtr),1))
		} else {
			rot.ranks[,1] <- rank(modtr)
			rot.ranks[,2] <- ngenes-rot.ranks[,1]+1
			rot.ranks[,3] <- rank(abs(modtr))
		}

		for(i in 1:nset)
		{
			rot.set.ranks <- colMeans(rot.ranks[iset[[i]],,drop=FALSE])
			p.value[i,] <- p.value[i,] + (rot.set.ranks >= obs.set.ranks[i,])
		}
	}	

	p.value <- (p.value+1)/(nrot+1)
	colnames(p.value)<-c("up","down","mixed")
	SetNames <- names(iset)
	if(is.null(SetNames))
		rownames(p.value) <- 1:nset
	else
		rownames(p.value) <- SetNames
	set.sizes <- unlist(lapply(iset,length))
	cbind(NGenes=set.sizes,p.value)
}

 