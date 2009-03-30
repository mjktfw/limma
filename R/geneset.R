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
# 24 Apr 2008. Last revised 9 Oct 2008.
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
	r1 <- mean(modt > sqrt(2))
	r2 <- mean(modt < -sqrt(2))
	out <- data.frame(Z=sqrt(statobs),Active=c(r1+r2,r1,r2,max(r1,r2)),P.Value=p)
	row.names(out) <- c("mixed","up","down","either")
	out
}

alias2Symbol <- function(alias,species="Hs",expand.symbols=FALSE)
#  Convert a set of alias names to official gene symbols
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  4 Sep 2008. Last revised 15 Jan 2009
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))
	DB <- paste("org",species,"eg","db",sep=".")
	ALIAS2EG <- paste("org",species,"egALIAS2EG",sep=".")
	SYMBOL <- paste("org",species,"egSYMBOL",sep=".")
	suppressPackageStartupMessages(require(DB,character.only=TRUE))
	if(expand.symbols)
	{
		alias <- intersect(alias,Rkeys(get(ALIAS2EG)))
		eg <- mappedLkeys(get(ALIAS2EG)[alias])
		mappedRkeys(get(SYMBOL)[eg])
	}
	else
	{
		isSymbol <- alias %in% Rkeys(get(SYMBOL)) 
		alias2 <- intersect(alias[!isSymbol],Rkeys(get(ALIAS2EG)))
		eg <- mappedLkeys(get(ALIAS2EG)[alias2])
		c(alias[isSymbol],mappedRkeys(get(SYMBOL)[eg]))

	}
}

barcodeplot <- function(selected,statistics,type="auto",...)
#	Barcode plot for gene set test
#	Gordon Smyth and Di Wu
#  20 October 2008. Last revised 21 Oct 2008.
{
	statistics <- as.numeric(statistics)
	isna <- is.na(statistics)
	if(any(isna)) {
		if(length(selected)==length(statistics)) {
			selected <- selected[!isna]
		} else {
			selected <- as.integer(selected)
			selected <- selected[!isna[selected]]
		}
		statistics <- statistics[!isna]
	}
	n <- length(statistics)
	type <- match.arg(type, c("t","f","auto"))
	allsamesign <- all(statistics >= 0) || all(statistics <= 0)
	if(type=="auto") {
		if(allsamesign)
			type <- "f"
		else
			type <- "t"
	}
	plot(1:n,xlim=c(0,n),ylim=c(0,1),type="n",axes=FALSE,xlab="",ylab="",...)
	npos <- sum(statistics > 1)
	nneg <- sum(statistics < -1)
	rect(npos+0.5,0,n-nneg+0.5,1,col="lightgray",border=NA)
	if(npos) rect(0.5,0,npos+0.5,1,col="pink",border=NA)
	if(nneg) rect(n-nneg+0.5,0,n+0.5,1,col="lightgreen",border=NA)
	r <- n+1-rank(statistics)[selected]
	lwd <- 50/length(r)
	lwd <- min(2,lwd)
	lwd <- max(0.1,lwd)
	segments(r,0,r,1,lwd=lwd)
#	rect(0.5,0,n+0.5,1,border="blue")
	if(type=="t") {
		mtext("Positive",side=2,line=-1,col="gray")
		mtext("Negative",side=4,line=-1,col="gray")
	} else {
		mtext("Largest",side=2,line=-1,col="gray")
		mtext("Smallest",side=4,line=-1,col="gray")
	}
	invisible()
}



##  ROMER.R

romer <- function(iset=NULL,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,nrot=10000)
# rotation-mean50-rank version of GSEA (gene set enrichment analysis) for linear models
# Gordon Smyth and Yifang Hu
# 30 March 2009.
{
	if(is.null(iset)) iset <- rep(TRUE,nrow(y))
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

# 	Fit model to all genes
	effects <- qr.qty(qr,t(y))
# 	Estimate global parameters s0 and d0
	s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
	sv <- squeezeVar(s2,df=d)
	d0 <- sv$df.prior
	s02 <- sv$var.prior
	sd.post <- sqrt(sv$var.post)

#	From here, all results are for set only
	Y <- effects[-(1:p0),,drop=FALSE]
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post
	
	s.r <- rank(modt)
	s.rank<-rep(0,nset)
	s.abs.rank<-rep(0,nset)
	s.either.rank<-rep(0,nset)
	s.rank.down<-rep(0,nset)

	modt.abs<-abs(modt)
	s.abs.r <-rank(modt.abs)

	m<-unlist(lapply(iset,length))
	m<-floor((m+1)/2)

	for(i in 1:nset)
	{	
		s.rank[i] <-.meanHalf(s.r[iset[[i]]],m[i])[1]
		s.abs.rank[i]<-.meanHalf(s.abs.r[iset[[i]]],m[i])[1]
		s.either.rank[i]<-.meanHalf(abs(s.r[iset[[i]]]-(ngenes+1)/2),m[i])[1]
		s.rank.down[i]<-.meanHalf(s.abs.r[iset[[i]]],m[i])[2]
	}

	p.value<-matrix(rep(0,nset*4),nrow=nset,ncol=4)

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
			if(.meanHalf(s.abs.r.2[iset[[j]]],m[j])[1]>=s.abs.rank[j]) p.value[j,1]<-p.value[j,1]+1
			if(.meanHalf(s.r.2[iset[[j]]],m[j])[1]>=s.rank[j]) p.value[j,2]<-p.value[j,2]+1
			if(.meanHalf(s.r.2[iset[[j]]],m[j])[2]>=s.rank.down[j]) p.value[j,3]<-p.value[j,3]+1
			if(.meanHalf(abs(s.r.2[iset[[j]]]-(ngenes+1)/2),m[j])[1]>=s.either.rank[j]) p.value[j,4]<-p.value[j,4]+1
		}
	}	

	p.value<-p.value/nrot
	colnames(p.value)<-c("mixed","up","down","either")
	rownames(p.value)<-names(iset)
	p.value
}

## Return means of top half and bottom half of the ranks for romer
.meanHalf<- function(x,n)
{
	l<-length(x)
	a<-sort(x,partial=n)
	top<-mean(a[1:n])
	if((l%%2)==0) bottom<-mean(a[(n+1):l])
	if((l%%2)!=0) bottom<-mean(a[n:l])
	c(top,bottom)
}


## MAKEIDX.R
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