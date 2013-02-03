##  ROMER.R

symbols2indices <- function(gene.sets, symbols, remove.empty=TRUE)
# Make a list of gene sets symbols into a list of gene sets indices used to create input for romer
# Gordon Smyth and Yifang Hu
# 25 March 2009.  Last modified 21 March 2011.
{
	gene.sets <- as.list(gene.sets)
	index <- lapply(gene.sets, function(x) which(symbols %in% x))
	if(remove.empty)
		for (i in length(index):1) if(!length(index[[i]])) index[[i]] <- NULL
	index
}

romer <- function(index,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation=NULL,set.statistic="mean",nrot=9999)
# rotation mean-rank version of GSEA (gene set enrichment analysis) for linear models
# Gordon Smyth and Yifang Hu
# 27 March 2009.  Last modified 3 June 2010.
{
	set.statistic <- match.arg(set.statistic,c("mean","floormean","mean50"))
	if(set.statistic=="mean50") {
		return(.romer.mean50(index=index,y=y,design=design,contrast=contrast,array.weights=array.weights,block=block,correlation=correlation,nrot=nrot))
	} else {
		return(.romer.mean.floormean(index=index,y=y,design=design,contrast=contrast,array.weights=array.weights,block=block,correlation=correlation,floor=(set.statistic=="floormean"),nrot=nrot))
	}
}

.romer.mean50 <- function(index,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,nrot=9999)
# rotation-mean50-rank version of GSEA (gene set enrichment analysis) for linear models
# Gordon Smyth and Yifang Hu
# 27 March 2009.  Last modified 15 Sep 2009.
{
#	Check input arguments
	if(!is.list(index)) index <- list(set=index)
	nset<-length(index)
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

	modt.abs<-abs(modt)
	s.abs.r <-rank(modt.abs)

	m<-unlist(lapply(index,length))
	m<-floor((m+1)/2)

	for(i in 1:nset)
	{	
		mh<-.meanHalf(s.r[index[[i]]],m[i])
		s.rank.up[i] <-mh[2]	
		s.rank.down[i]<-mh[1]
  		s.rank.mixed[i]<-.meanHalf(s.abs.r[index[[i]]],m[i])[2]
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
	p.value <- matrix(rep(0,nset*3),nrow=nset,ncol=3)
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
			mh.2<-.meanHalf(s.r.2[index[[j]]],m[j])
			
			s.rank.up.2 <-mh.2[2]	
			s.rank.down.2 <-mh.2[1]
  			s.rank.mixed.2 <-.meanHalf(s.abs.r.2[index[[j]]],m[j])[2]
		
			if(s.rank.mixed.2>=s.rank.mixed[j]) p.value[j,1]<-p.value[j,1]+1
			if(s.rank.up.2>=s.rank.up[j]) p.value[j,2]<-p.value[j,2]+1
			if(s.rank.down.2<=s.rank.down[j]) p.value[j,3]<-p.value[j,3]+1
		}
	}	

	p.value <- (p.value+1)/(nrot+1)
	colnames(p.value)<-c("Mixed","Up","Down")
	SetNames <- names(index)
	if(is.null(SetNames))
		rownames(p.value) <- 1:nset
	else
		rownames(p.value) <- SetNames
	len.index<-as.numeric(lapply(index,length))
	cbind(NGenes=len.index,p.value)
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

.romer.mean.floormean <- function(index,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,floor=FALSE,nrot=9999)
# rotation mean-rank version of GSEA (gene set enrichment analysis) for linear models
# Gordon Smyth and Yifang Hu
# 27 March 2009.  Last modified 5 Oct 2009.
{
#	Check input arguments
	if(!is.list(index)) index <- list(set=index)
	nset<-length(index)
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
	for(i in 1:nset) obs.set.ranks[i,] <- colMeans(obs.ranks[index[[i]],,drop=FALSE])

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

		for(k in 1:nset)
		{
			rot.set.ranks <- colMeans(rot.ranks[index[[k]],,drop=FALSE])
			p.value[k,] <- p.value[k,] + (rot.set.ranks >= obs.set.ranks[k,])
		}
	}	

	p.value <- (p.value+1)/(nrot+1)
	colnames(p.value)<-c("Up","Down","Mixed")
	SetNames <- names(index)
	if(is.null(SetNames))
		rownames(p.value) <- 1:nset
	else
		rownames(p.value) <- SetNames
	set.sizes <- unlist(lapply(index,length))
	cbind(NGenes=set.sizes,p.value)
}

topRomer<-function(x,n=10,alternative="up")
# extracts a number of top gene sets results from the romer output.
# Gordon Smyth and Yifang Hu.
# 22 Mar 2010.  Last modified 5 Sep 2011.
{
	n <- min(n,nrow(x))
	alternative <- match.arg(tolower(alternative),c("up","down","mixed"))
	alternative <- switch(alternative,"up"="Up","down"="Down","mixed"="Mixed")
	o <- switch(alternative,
		"Up"=order(x[,"Up"],x[,"Mixed"],-x[,"NGenes"]),
		"Down"=order(x[,"Down"],x[,"Mixed"],-x[,"NGenes"]),
		"Mixed"=order(x[,"Mixed"],x[,"Up"],x[,"Down"],-x[,"NGenes"])
	)
	x[o,][1:n,,drop=FALSE]
}

