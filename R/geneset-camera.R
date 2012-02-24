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


camera <- function(indices,y,design,contrast=ncol(design),statistic="modt",use.ranks=FALSE,trend.var=FALSE)
#	Competitive gene set test allowing for correlation between genes
#	Gordon Smyth and Di Wu
#  Created 1 Feb 2012.  Last modified 24 Feb 2012
{
	if(!is.list(indices)) indices <- list(set1=indices)
	y <- as.matrix(y)
	statistic <- match.arg(statistic,choices=c("modt","logFC"))

	if(is.character(contrast)) {
		contrast <- which(contrast==colnames(design))
		if(length(contrast)==0) stop("coef ",contrast," not found")
	}
#	Reform design matrix so that contrast of interest is last column
	coef <- p <- ncol(design)
	if(length(contrast)==1) {
		if(contrast<p) design <- cbind(design[,-contrast,drop=FALSE],design[,contrast])
	} else {
		qr <- qr(contrast)
		Q <- qr.Q(qr,complete=TRUE)
		Q <- cbind(Q[,-1],sign(qr$qr[1,1])*Q[,1])
		design <- design %*% Q
	}

	G <- nrow(y)
	fit <- lm.fit(design,t(y))
	df.camera <- min(fit$df.residual,G-2)
	U <- fit$effects[-(1:fit$qr$rank),,drop=FALSE]
	sigma2 <- colMeans(U^2)
	U <- t(U) / sqrt(sigma2)
	if(statistic=="logFC") {
		Stat <- fit$coefficients[coef,]
	} else {
		if(trend.var) A <- rowMeans(y) else A <- NULL
		sv <- squeezeVar(sigma2,df=fit$df.residual,covariate=A)
		modt <- sign(fit$qr$qr[coef,coef]) * fit$effects[coef,] / sqrt(sv$var.post)
		df.total <- pmin(fit$df.residual+sv$df.prior, G*fit$df.residual)
		Stat <- zscoreT(modt, df=df.total)
	}
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
