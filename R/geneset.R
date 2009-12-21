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

wilcoxGST <- function(selected,statistics,alternative="mixed")
#	Mean-rank gene set test using Wilcox test.
#	Gordon Smyth
#	27 July 2009.  Last modified 28 July 2009.

geneSetTest(selected=selected,statistics=statistics,alternative=alternative,type="t",ranks.only=TRUE)


roast <- function(iset=NULL,y,design,contrast=ncol(design),gene.weights=NULL,array.weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,nrot=999)
# rotation gene set testing for linear models
# Gordon Smyth and Di Wu
# 24 Apr 2008. Last revised 15 May 2009.
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

#	Reform design matrix so that contrast of interest is last column
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
	p <- (rowSums(t(stati) >= statobs)+1)/(nrot+1)

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

barcodeplot <- function(selected,statistics,labels=c("Up","Down"),...)
#	Barcode plot for gene set test
#	Gordon Smyth and Di Wu
#  20 October 2008. Last revised 30 Sep 2009.
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
	mtext(labels[1],side=2,line=-1,col="gray")
	mtext(labels[2],side=4,line=-1,col="gray")
	invisible()
}

barcodeplot2 <- function (selected,statistics,selected2=NULL,labels=c("Up","Down"), ...) 
#	Barcode plot for gene set test. 
#  Plots both up and down sets in one figure.
#	Gordon Smyth and Di Wu
#  20 October 2008. Last revised 30 Sep 2009.
{
	statistics <- as.numeric(statistics)
	isna <- is.na(statistics)
	if (any(isna)) {
		if (length(selected) == length(statistics)) {
			selected <- selected[!isna]
		}
		else {
			selected <- as.integer(selected)
			selected <- selected[!isna[selected]]
		}
		statistics <- statistics[!isna]
	}
	n <- length(statistics)

	plot(1:n,xlim=c(0,n),ylim=c(-1,1),type="n",axes=FALSE,xlab="",ylab="",...)
	npos <- sum(statistics > 1)
	nneg <- sum(statistics < -1)
	rect(npos+0.5,-0.5,n-nneg+0.5,0.5,col="lightgray",border=NA)
	if(npos) rect(0.5,-0.5,npos+0.5,0.5,col="pink",border=NA)
	if(nneg) rect(n-nneg+0.5,-0.5,n+0.5,0.5,col="lightblue",border=NA)

	rankstat <- rank(statistics)
	r <- n+1-rankstat[selected]
	lwd <- 50/length(r)
	lwd <- min(2,lwd)
	lwd <- max(0.1,lwd)
	segments(r,0,r,1,lwd=lwd,col="red")

	if(!is.null(selected2)) {
		r2 <- n+1-rankstat[selected2]
		lwd2 <- 50/length(r2)
		lwd2 <- min(2,lwd2)
		lwd2 <- max(0.1,lwd2)
		segments(r2,-1,r2,0,lwd=lwd2,col="blue")
	}

	mtext(labels[1],side=2,line=-0.5,col="black")
	mtext(labels[2],side=4,line=-1,col="black")
	invisible()
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


alias2SymbolTable <- function(alias,species="Hs")
#  Convert a set of alias names to official gene symbols of the same length
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  3 Sep 2009.  Last modified 17 Dec 2009.
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))
	DB <- paste("org",species,"eg","db",sep=".")
	ALIAS2EG <- paste("org",species,"egALIAS2EG",sep=".")
	SYMBOL <- paste("org",species,"egSYMBOL",sep=".")
	suppressPackageStartupMessages(require(DB,character.only=TRUE))
	
	isSymbol <- alias %in% Rkeys(get(SYMBOL)) 
	Symbol<-rep.int(NA,length(alias))
	Symbol[isSymbol]<-alias[isSymbol]
				
	isalias<-(alias[!isSymbol]) %in% (Rkeys(get(ALIAS2EG)))
	alias2<-(alias[!isSymbol])[isalias]

	aliasTbl<-toTable(get(ALIAS2EG)[alias2])
	hits<-names(table(aliasTbl$alias_symbol))[as.numeric(table(aliasTbl$alias_symbol))>1]
	if(length(hits)>0) warning("Multiple Hits for ", hits)
	
	aliasTbl.o<-aliasTbl[match(alias2,aliasTbl$alias_symbol),]
	symb<-toTable(get(SYMBOL)[aliasTbl.o$gene_id])
	Symbol[!isSymbol][isalias]<-symb[match(aliasTbl.o$gene_id,symb$gene_id),]$symbol
	Symbol	
}

##  ROMER2.R

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


##  ROMER.R

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
