##  CAMERA.R

interGeneCorrelation <- function(y, design)
#	Estimate variance-inflation factor for means of correlated genes
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 11 Feb 2012
{
	m <- nrow(y)
	qrdesign <- qr(design)
	y <- qr.qty(qrdesign, t(y))[-(1:qrdesign$rank),]
#	Gives same result as the following
#	ny <- t(y) / sqrt(colSums(y^2))
#	cormatrix <- tcrossprod(ny)
#	correlation <- mean(cormatrix[lower.tri(cormatrix)])
#	1+correlation*(n-1)
	y <- t(y) / sqrt(colMeans(y^2))
	vif <- m * mean(colMeans(y)^2)
	correlation <- (vif-1)/(m-1)
	list(vif=vif,correlation=correlation)
}


camera <- function(indices,y,design,contrast=ncol(design),weights=NULL,use.ranks=FALSE,trend.var=FALSE)
#	Competitive gene set test allowing for correlation between genes
#	Gordon Smyth and Di Wu
#  Created 2007.  Last modified 26 Feb 2012
{
#	Check arguments
	if(!is.list(indices)) indices <- list(set1=indices)
	if(is(y,"EList") || is(y,"MAList")) {
		if(is.null(design)) design <- as.matrix(y$design)
		if(is.null(weights)) weights <- as.matrix(y$weights)
	}
	y <- as.matrix(y)

	G <- nrow(y)
	n <- ncol(y)
	p <- ncol(design)
	df.residual <- n-p
	df.camera <- min(df.residual,G-2)

#	Reform design matrix so that contrast of interest is last column
	if(is.character(contrast)) {
		contrast <- which(contrast==colnames(design))
		if(length(contrast)==0) stop("coef ",contrast," not found")
	}
	if(length(contrast)==1) {
		j <- c((1:p)[-contrast], contrast)
		if(contrast<p) design <- design[,j]
	} else {
		QR <- qr(contrast)
		design <- t(qr.qty(QR,t(design)))
		if(sign(QR$qr[1,1]<0)) design[,1] <- -design[,1]
		design <- design[,c(2:p,1)]
	}

#	Compute effects matrix
	if(is.null(weights)) {
		QR <- qr(design)
		if(QR$rank<p) stop("design matrix is not of full rank")
		effects <- qr.qty(QR,t(y))
		unscaledt <- effects[p,]
		if(QR$qr[p,p]<0) unscaledt <- -unscaledt
	} else {
		if(any(weights<=0)) stop("only positive weights permitted")
		effects <- matrix(0,n,G)
		unscaledt <- rep(0,n)
		sw <- sqrt(weights)
		yw <- y*sw
		for (g in 1:G) {
			xw <- design*sw[g,]
			QR <- qr(xw)
			if(QR$rank<p) stop("weighted design matrix not of full rank for gene ",g)
			effects[,g] <- qr.qty(QR,yw[g,])
			unscaledt[g] <- effects[p,g]
			if(QR$qr[p,p]<0) unscaledt[g] <- -unscaledt[g]
		}
	}

#	Standardized residuals
	U <- effects[-(1:p),,drop=FALSE]
	sigma2 <- colMeans(U^2)
	U <- t(U) / sqrt(sigma2)

#	Moderated t
	if(trend.var) A <- rowMeans(y) else A <- NULL
	sv <- squeezeVar(sigma2,df=df.residual,covariate=A)
	modt <- unscaledt / sqrt(sv$var.post)
	df.total <- min(df.residual+sv$df.prior, G*df.residual)
	Stat <- zscoreT(modt, df=df.total)

#	Global statistics
	meanStat <- mean(Stat)
	varStat <- var(Stat)

	nsets <- length(indices)
	tab <- matrix(0,nsets,5)
	rownames(tab) <- names(indices)
	colnames(tab) <- c("NGenes","Correlation","Down","Up","TwoSided")
	for (i in 1:nsets) {
		index <- indices[[i]]
		StatInSet <- Stat[index]
		m <- length(StatInSet)
		m2 <- G-m
		if(m>1) {
			Uset <- U[index,,drop=FALSE]
			vif <- m * mean(colMeans(Uset)^2)
			correlation <- (vif-1)/(m-1)
		} else {
			vif <- 1
			correlation <- NA
		}

		if(use.ranks) {
			tab[i,3:4] <- rankSumTestWithCorrelation(index,statistics=Stat,correlation=correlation,df=df.camera)
		} else {	
			meanStatInSet <- mean(StatInSet)
			delta <- G/m2*(meanStatInSet-meanStat)
			varStatPooled <- ( (G-1)*varStat - delta^2*m*m2/G ) / (G-2)
			two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
			tab[i,3] <- pt(two.sample.t,df=df.camera)
			tab[i,4] <- pt(two.sample.t,df=df.camera,lower.tail=FALSE)
		}
		tab[i,1] <- m
		tab[i,2] <- correlation
	}
	tab[,5] <- 2*pmin(tab[,3],tab[,4])
	tab
}
