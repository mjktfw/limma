fitFDistRobustly <- function(x,df1,covariate=NULL,winsor.tail.p=c(0.05,0.1),trace=FALSE)
#	Robust estimation of the parameters of a scaled F-distribution,
#	given the first degrees of freedom, using first and second
#	moments of Winsorized z-values
#	Gordon Smyth and Belinda Phipson
#	8 Sept 2002.  Last revised 14 January 2015.
{
#	Check x
	n <- length(x)

#	Eliminate case of no useful data
	if(n==0) 
		if(is.null(covariate))
			return(list(scale=NaN,df2=NaN,df2.robust=NaN))
		else
			return(list(scale=numeric(0),df2=numeric(0),df2.robust=NaN))

#	Check df1
	if(length(df1)>1) if(length(df1) != n) stop("x and df1 are different lengths")

#	Check covariate
	if(!is.null(covariate)) {
		if(length(covariate) != n) stop("x and covariate are different lengths")
		if(!all(is.finite(covariate))) stop("covariate contains NA or infinite values")
	}
	
#	Treat zero df1 values as non-informative cases
#	Similarly for missing values or x or df1
	ok <- !is.na(x) & is.finite(df1) & (df1 > 1e-6)
	notallok <- !all(ok)
	if(notallok) {
		df2.shrunk <- x
		x <- x[ok]
		if(length(df1)>1) df1 <- df1[ok]
		if(!is.null(covariate)) {
			covariate2 <- covariate[!ok]
			covariate <- covariate[ok]
		}
		fit <- Recall(x=x,df1=df1,covariate=covariate,winsor.tail.p=winsor.tail.p,trace=trace)
		df2.shrunk[ok] <- fit$df2.shrunk
		df2.shrunk[!ok] <- fit$df2
		if(is.null(covariate))
			scale <- fit$scale
		else {
			scale <- x
			scale[ok] <- fit$scale
			scale[!ok] <- exp(approx(covariate,log(fit$scale),xout=covariate2,rule=2)$y)
		}
		return(list(scale=scale,df2=fit$df2,df2.shrunk=df2.shrunk))
	}

#	Avoid zero or negative x values
	m <- median(x)
	if(m<=0)	stop("x values are mostly <= 0")
	i <- (x < m*1e-12)
	if(any(i)) {
		warning("small x values have been offset away from zero")
		x[i] <- m*1e-12
	}

#	Store non-robust estimates
	NonRobust <- fitFDist(x=x,df1=df1,covariate=covariate)

#	Check winsor.tail.p
	prob <- winsor.tail.p <- rep(winsor.tail.p,length=2)
	prob[2] <- 1-winsor.tail.p[2]
	if(all(winsor.tail.p < 1/n)) {
		NonRobust$df2.shrunk <- rep.int(NonRobust$df2,n)
		return(NonRobust)
	}

#	Transform x to constant df1
	if(length(df1)>1) {
		df1max <- max(df1)
		i <- (df1 < (df1max - 1e-14))
		if(any(i)) {
			if(is.null(covariate)) s <- NonRobust$scale else s <- NonRobust$scale[i]
			f <- x[i]/s
			df2 <- NonRobust$df2
			pupper <- pf(f,df1=df1[i],df2=df2,lower.tail=FALSE,log.p=TRUE)
			plower <- pf(f,df1=df1[i],df2=df2,lower.tail=TRUE,log.p=TRUE)
			up <- pupper<plower
			if(any(up)) f[up] <- qf(pupper[up],df1=df1max,df2=df2,lower.tail=FALSE,log.p=TRUE)
			if(any(!up)) f[!up] <- qf(plower[!up],df1=df1max,df2=df2,lower.tail=TRUE,log.p=TRUE)
			x[i] <- f*s
			df1 <- df1max
		} else {
			df1 <- df1[1]
		}
	}

#	Better to work with log(F)
	z <- log(x)

#	Demean or Detrend
	if(is.null(covariate)) {
		ztrend <- mean(z,trim=winsor.tail.p[2])
		zresid <- z-ztrend
	} else {
		lo <- loessFit(z,covariate,span=0.4)
		ztrend <- lo$fitted
		zresid <- lo$residual
	}

#	Moments of Winsorized residuals
	zrq <- quantile(zresid,prob=prob)
	zwins <- pmin(pmax(zresid,zrq[1]),zrq[2])
	zwmean <- mean(zwins)
	zwvar <- mean((zwins-zwmean)^2)*n/(n-1)
	if(trace) cat("Variance of Winsorized Fisher-z",zwvar,"\n")

#	Theoretical Winsorized moments
	if(!requireNamespace("statmod",quietly=TRUE)) stop("statmod package required but is not available")
	g <- statmod::gauss.quad.prob(128,dist="uniform")
	linkfun <- function(x) x/(1+x)
	linkinv <- function(x) x/(1-x)
	winsorizedMoments <- function(df1=df1,df2=df2,winsor.tail.p=winsor.tail.p) {
		fq <- qf(p=c(winsor.tail.p[1],1-winsor.tail.p[2]),df1=df1,df2=df2)
		zq <- log(fq)
		q <- linkfun(fq)
		nodes <- q[1] + (q[2]-q[1]) * g$nodes
		fnodes <- linkinv(nodes)
		znodes <- log(fnodes)
		f <- df(fnodes,df1=df1,df2=df2)/(1-nodes)^2
		q21 <- q[2]-q[1]
		m <- q21*sum(g$weights*f*znodes) + sum(zq * winsor.tail.p)
		v <- q21*sum(g$weights*f*(znodes-m)^2) + sum((zq-m)^2 * winsor.tail.p)
		list(mean=m,var=v)
	}

#	Try df2==Inf
	mom <- winsorizedMoments(df1=df1,df2=Inf,winsor.tail.p=winsor.tail.p)
	funvalInf <- log(zwvar/mom$var)
	if(funvalInf <= 0) {
		df2 <- Inf
#		Correct trend for bias
		ztrendcorrected <- ztrend+zwmean-mom$mean
		s20 <- exp(ztrendcorrected)
#		Posterior df for outliers
		Fstat <- exp(z-ztrendcorrected)
		TailP <- pchisq(Fstat*df1,df=df1,lower.tail=FALSE)
		r <- rank(Fstat)
		EmpiricalTailProb <- (n-r+0.5)/n
		ProbNotOutlier <- pmin(TailP/EmpiricalTailProb,1)
		df.pooled <- n*df1
		df2.shrunk <- rep.int(df2,n)
		O <- ProbNotOutlier < 1
		if(any(O)) {
			df2.shrunk[O] <- ProbNotOutlier[O]*df.pooled
			o <- order(TailP)
			df2.shrunk[o] <- cummax(df2.shrunk[o])
		}
		return(list(scale=s20,df2=df2,df2.shrunk=df2.shrunk))
	}

#	Estimate df2 by matching variance of zwins
#	Use beta distribution Gaussian quadrature to find mean and variance
#	of Winsorized F-distribution
	fun <- function(x) {
		df2 <- linkinv(x)
		mom <- winsorizedMoments(df1=df1,df2=df2,winsor.tail.p=winsor.tail.p)
		if(trace) cat("df2=",df2,", Working Var=",mom$var,"\n")
		log(zwvar/mom$var)
	}

#	Use non-robust estimate as lower bound for df2
	if(NonRobust$df2==Inf) {
		NonRobust$df2.shrunk <- rep.int(NonRobust$df2,n)
		return(NonRobust)
	}
	rbx <- linkfun(NonRobust$df2)
	funvalLow <- fun(rbx)
	if(funvalLow >= 0) {
		df2 <- NonRobust$df2
	} else {
		u <- uniroot(fun,interval=c(rbx,1),tol=1e-8,f.lower=funvalLow,f.upper=funvalInf)
		df2 <- linkinv(u$root)
	}

#	Correct ztrend for bias
	mom <- winsorizedMoments(df1=df1,df2=df2,winsor.tail.p=winsor.tail.p)
	ztrendcorrected <- ztrend+zwmean-mom$mean
	s20 <- exp(ztrendcorrected)

#	Posterior df for outliers
	zresid <- z-ztrendcorrected
	Fstat <- exp(zresid)
	TailP <- pf(Fstat,df1=df1,df2=df2,lower.tail=FALSE)
	r <- rank(Fstat)
	EmpiricalTailProb <- (n-r+0.5)/n
	ProbNotOutlier <- pmin(TailP/EmpiricalTailProb,1)
	df2.shrunk <- rep.int(df2,n)
	if(any(ProbNotOutlier<1)) {
		ProbOutlier <- 1-ProbNotOutlier
		VarOutlier <- max(zresid)^2
		VarOutlier <- VarOutlier-trigamma(df1/2)
		if(trace) cat("VarOutlier",VarOutlier,"\n")
		if(VarOutlier > 0) {
			df2.outlier <- 2*trigammaInverse(VarOutlier)
			if(trace) cat("df2.outlier",df2.outlier,"\n")
			if(df2.outlier < df2) {
				df2.shrunk <- ProbNotOutlier*df2+ProbOutlier*df2.outlier
				o <- order(TailP)
				df2.ordered <- df2.shrunk[o]
				df2.ordered[1] <- min(df2.ordered[1],NonRobust$df2)
				m <- cumsum(df2.ordered)
				m <- m/(1:n)
				imin <- which.min(m)
				df2.ordered[1:imin] <- m[imin]
				df2.shrunk[o] <- cummax(df2.ordered)
			}
		}
	}

	list(scale=s20,df2=df2,df2.shrunk=df2.shrunk)
}
