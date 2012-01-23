##  CAMERA.R

interGeneCorrelation <- function(y, design)
#	Estimate average correlation between genes
#	Returns the variance-inflation factor 1+correlation*(n-1)
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 19 Jan 2012
{
	qrdesign <- qr(design)
	y <- qr.qty(qrdesign, t(y))[-(1:qrdesign$rank),]
#	Gives same result as the following
#	ny <- t(y) / sqrt(colSums(y^2))
#	cormatrix <- tcrossprod(ny)
#	correlation <- mean(cormatrix[lower.tri(cormatrix)])
	y <- t(y) / sqrt(colMeans(y^2))
	nrow(y) * mean(colMeans(y)^2)
}


rankSumTestwithCorrelation <- function(index,statistics,correlation=NULL,vif=NULL,df=Inf)
#	Rank sum test as for two-sample Wilcoxon-Mann-Whitney test,
#	but allowing for inter-gene correlation in test set
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 20 Jan 2012.
{
	n <- length(statistics)
	r <- rank(statistics)
	r1 <- r[index]
	n1 <- length(r1)
	n2 <- n-n1
	U <- n1*n2 + n1*(n1+1)/2 - sum(r1)
	mu <- n1*n2/2
	TIES <- (length(r) != length(unique(r)))
	if(TIES) {
		NTIES <- table(r)
		sumt <- sum(NTIES*(NTIES+1)*(NTIES-1))
		sigma2 <- n1*n2/n/(n-1)*(n*(n+1)*(n-1)-sumt)/12
	} else {
		sigma2 <- n1*n2*(n+1)/12
	}
#	Variance inflation
	if(is.null(vif) && !is.null(correlation)) vif <- 1+(n1-1)*correlation*(n-n1)/n
	if(!is.null(vif)) sigma2 <- sigma2 * vif
	zlowertail <- (U+0.5-mu)/sqrt(sigma2)
	zuppertail <- (U-0.5-mu)/sqrt(sigma2)

#	Lower and upper tails are reversed on output
#	because R's ranks are the reverse of Mann-Whitney's ranks
	pvalues <- c(lower.tail=pt(zuppertail,df=df,lower.tail=FALSE), upper.tail=pt(zlowertail,df=df))
	pvalues	
}


camera <- function(index,y,design,coef=ncol(design),statistic="modt",use.ranks=TRUE,correlation=NULL,vif=NULL,df=Inf)
#	Competitive gene set test allowing for correlation between genes
#	Di Wu and Gordon Smyth
#  Created 2007.  Last modified 21 Jan 2012
{
	statistic <- match.arg(statistic,c("modt","logFC"))
	fit <- lmFit(y,design)
	if(statistic=="logFC") {
		Stat <- fit$coef[,coef]
	} else {
		fit <- eBayes(fit)
		Stat <- zscoreT(fit$t[,coef],df=fit$df.total)
	}
	G <- length(Stat)
	StatInSet <- Stat[index]
	m <- length(StatInSet)

	if(is.null(vif) && !is.null(correlation)) vif <- 1+(m-1)*correlation
	if(is.null(vif)) {
		vif <- interGeneCorrelation(y[index,], design)
		df <- fit$df.residual[1]
	}
	vif <- vif*(G-m)/G

	if(use.ranks) {
		pvalues <- rankSumTestwithCorrelation(index,statistics=Stat,vif=vif,df=df)
	} else {
		mu <- mean(Stat)
		se2 <- var(Stat)/m
		Z <- (mean(StatInSet) - mu) / sqrt(se2*vif)
		pvalues <- c(lower.tail=pt(Z,df=df),upper.tail=pt(Z,df=df,lower=FALSE))
	}
	pvalues["two.tail"] <- 2*min(pvalues)
	list(p.values=pvalues,vif=vif,df=df)
}


cameraFaster <- function(index,y,design,coef=ncol(design),statistic="modt")
#	Competitive gene set test allowing for correlation between genes
#	Gordon Smyth and Di Wu
#  Created 2007.  Last modified 21 Jan 2012
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
	vif <- m * mean(colMeans(U)^2) * (G-m)/G

	if(statistic=="logFC") {
		Stat <- fit$coefficients[coef,]
	} else {
		sv <- squeezeVar(sigma2,df=fit$df.residual)
		modt <- fit$effects[coef,] / sqrt(sv$var.post)
		df.total <- pmin(fit$df.residual+sv$df.prior, G*fit$df.residual)
		Stat <- zscoreT(modt, df=df.total)
	}
	StatInSet <- Stat[index]

	df.camera <- fit$df.residual
	pvalues <- matrix(0,2,3)
	rownames(pvalues) <- c("Ranks","Parametric")
	colnames(pvalues) <- c("lower.tail","upper.tail","two.tail")

	pvalues["Ranks",1:2] <- rankSumTestwithCorrelation(index,statistics=Stat,vif=vif,df=df.camera)
	pvalues["Ranks","two.tail"] <- 2*min(pvalues[1,1:2])
	
	mu <- mean(Stat)
	se2 <- var(Stat)/m
	Z <- (mean(StatInSet) - mu) / sqrt(se2*vif)
	pvalues["Parametric","lower.tail"] <- pt(Z,df=df.camera)
	pvalues["Parametric","upper.tail"] <- pt(Z,df=df.camera,lower=FALSE)
	pvalues["Parametric","two.tail"] <- 2*min(pvalues[2,1:2])

	list(p.values=pvalues,vif=vif,df=df.camera)
}

