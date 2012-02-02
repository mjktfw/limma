##  CAMERA.R

camera <- function(index,y,design,coef=ncol(design),statistic="modt")
#	Competitive gene set test allowing for correlation between genes
#	Gordon Smyth and Di Wu
#  Created 2007.  Last modified 2 Feb 2012
{
	statistic <- match.arg(statistic,choices=c("modt","logFC"))
	fit <- lm.fit(design,t(y))
	U <- fit$effects[-(1:fit$qr$rank),,drop=FALSE]
	sigma2 <- colMeans(U^2)
	G <- nrow(y)

#	Estimate VIF
	U <- U[,index,drop=FALSE]
	U <- t(U) / sqrt(sigma2[index])
	m <- nrow(U)
	vif <- m * mean(colMeans(U)^2)

	if(statistic=="logFC") {
		Stat <- fit$coefficients[coef,]
	} else {
		sv <- squeezeVar(sigma2,df=fit$df.residual)
		modt <- sign(fit$qr$qr[coef,coef]) * fit$effects[coef,] / sqrt(sv$var.post)
		df.total <- pmin(fit$df.residual+sv$df.prior, G*fit$df.residual)
		Stat <- zscoreT(modt, df=df.total)
	}
	StatInSet <- Stat[index]

	df.camera <- min(fit$df.residual,G-2)
	pvalues <- matrix(0,2,3)
	rownames(pvalues) <- c("Ranks","Parametric")
	colnames(pvalues) <- c("lower.tail","upper.tail","two.tail")
	pvalues <- pvalues

	correlation <- (vif-1)/(m-1)
	pvalues["Ranks",1:2] <- rankSumTestwithCorrelation(index,statistics=Stat,correlation=correlation,df=df.camera)
	pvalues["Ranks","two.tail"] <- 2*min(pvalues["Ranks",1:2])
	
	mu <- mean(Stat)
	s2 <- var(Stat)
	meanStat <- mean(StatInSet)
	m2 <- G-m
	delta <- G/m2*(meanStat-mu)
	s2pool <- ( (G-1)*s2 - delta^2*m*m2/G ) / (G-2)
	two.sample.t <- delta / sqrt( s2pool * (vif/m + 1/m2) )
	pvalues["Parametric","lower.tail"] <- pt(two.sample.t,df=df.camera)
	pvalues["Parametric","upper.tail"] <- pt(two.sample.t,df=df.camera,lower=FALSE)
	pvalues["Parametric","two.tail"] <- 2*min(pvalues[2,1:2])

	list(p.values=pvalues,vif=vif,correlation=correlation,df=df.camera)
}


interGeneVIF <- function(y, design)
#	Estimate variance-inflation factor for means of correlated genes
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 30 Jan 2012
{
	qrdesign <- qr(design)
	y <- qr.qty(qrdesign, t(y))[-(1:qrdesign$rank),]
#	Gives same result as the following
#	ny <- t(y) / sqrt(colSums(y^2))
#	cormatrix <- tcrossprod(ny)
#	correlation <- mean(cormatrix[lower.tri(cormatrix)])
#	1+correlation*(n-1)
	y <- t(y) / sqrt(colMeans(y^2))
	nrow(y) * mean(colMeans(y)^2)
}


rankSumTestwithCorrelation <- function(index,statistics,correlation=0,df=Inf)
#	Rank sum test as for two-sample Wilcoxon-Mann-Whitney test,
#	but allowing for inter-gene correlation in test set
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 29 Jan 2012.
{
	n <- length(statistics)
	r <- rank(statistics)
	r1 <- r[index]
	n1 <- length(r1)
	n2 <- n-n1
	U <- n1*n2 + n1*(n1+1)/2 - sum(r1)
	mu <- n1*n2/2

	if(correlation==0 ) {
		sigma2 <- n1*n2*(n+1)/12
	} else {
		sigma2 <- asin(1)*n1*n2 + asin(0.5)*n1*n2*(n2-1) + asin(correlation/2)*n1*(n1-1)*n2*(n2-1) + asin((correlation+1)/2)*n1*(n1-1)*n2
		sigma2 <- sigma2/2/pi
	}

	TIES <- (length(r) != length(unique(r)))
	if(TIES) {
		NTIES <- table(r)
		adjustment <- sum(NTIES*(NTIES+1)*(NTIES-1)) / (n*(n+1)*(n-1))
		sigma2 <- sigma2 * (1 - adjustment)
	}
	zlowertail <- (U+0.5-mu)/sqrt(sigma2)
	zuppertail <- (U-0.5-mu)/sqrt(sigma2)

#	Lower and upper tails are reversed on output
#	because R's ranks are the reverse of Mann-Whitney's ranks
	pvalues <- c(lower.tail=pt(zuppertail,df=df,lower.tail=FALSE), upper.tail=pt(zlowertail,df=df))
	pvalues	
}


mcamera <- function(indices,y,design,coef=ncol(design),statistic="modt")
#	Competitive gene set test allowing for correlation between genes
#	Gordon Smyth and Di Wu
#  Created 1 Feb 2012.  Last modified 2 Feb 2012
{
	if(!is.list(indices)) indices <- list(set1=indices)
	statistic <- match.arg(statistic,choices=c("modt","logFC"))
	G <- nrow(y)
	fit <- lm.fit(design,t(y))
	df.camera <- min(fit$df.residual,G-2)
	U <- fit$effects[-(1:fit$qr$rank),,drop=FALSE]
	sigma2 <- colMeans(U^2)
	U <- t(U) / sqrt(sigma2)
	if(statistic=="logFC") {
		Stat <- fit$coefficients[coef,]
	} else {
		sv <- squeezeVar(sigma2,df=fit$df.residual)
		modt <- sign(fit$qr$qr[coef,coef]) * fit$effects[coef,] / sqrt(sv$var.post)
		df.total <- pmin(fit$df.residual+sv$df.prior, G*fit$df.residual)
		Stat <- zscoreT(modt, df=df.total)
	}
	meanStat <- mean(Stat)
	varStat <- var(Stat)

	nsets <- length(indices)
	tab <- matrix(0,nsets,5)
	rownames(tab) <- names(indices)
	colnames(tab) <- c("NGenes","Correlation","P.Lower","P.Upper","P.TwoTail")
	for (i in 1:nsets) {
		index <- indices[[i]]
		Uset <- U[index,,drop=FALSE]
		m <- nrow(Uset)
		m2 <- G-m
		if(m>1) {
			vif <- m * mean(colMeans(Uset)^2)
			correlation <- (vif-1)/(m-1)
		} else {
			vif <- 1
			correlation <- NA
		}
		StatInSet <- Stat[index]
		meanStatInSet <- mean(StatInSet)
		delta <- G/m2*(meanStatInSet-meanStat)
		varStatPooled <- ( (G-1)*varStat - delta^2*m*m2/G ) / (G-2)
		two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
		tab[i,1] <- m
		tab[i,2] <- correlation
		tab[i,3] <- pt(two.sample.t,df=df.camera)
		tab[i,4] <- pt(two.sample.t,df=df.camera,lower=FALSE)
	}
	tab[,5] <- 2*pmin(tab[,3],tab[,4])
	tab
}
