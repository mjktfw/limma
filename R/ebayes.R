#  DIFFERENTIAL EXPRESSION

eBayes <- function(fit,proportion=0.01,std.coef=NULL) {
#	Empirical Bayes statistics to select differentially expressed genes
#	Object orientated version
#	Gordon Smyth
#	4 August 2003.

	eb <- ebayes(fit=fit,proportion=proportion,std.coef=std.coef)
	fit$df.prior <- eb$df.prior
	fit$s2.prior <- eb$s2.prior
	fit$var.prior <- eb$var.prior
	fit$proportion <- proportion
	fit$s2.post <- eb$s2.post
	fit$t <- eb$t
	fit$p.value <- eb$p.value
	fit$lods <- eb$lods
	fit
}

ebayes <- function(fit,proportion=0.01,std.coef=NULL) {
#	Empirical Bayes statistics to select differentially expressed genes
#	Gordon Smyth
#	8 Sept 2002.  Last revised 14 August 2003.

	coefficients <- fit$coefficients
	stdev.unscaled <- fit$stdev.unscaled
	sigma <- fit$sigma
	df.residual <- fit$df.residual
	if(is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) stop("No data, or argument is not a valid lmFit object")
	if(all(df.residual==0)) stop("No residual degrees of freedom in linear model fits")
	if(all(!is.finite(sigma))) stop("No finite residual standard deviations")

#	Moderated t-statistic
	out <- fitFDist(sigma^2,df1=df.residual)
	out$s2.prior <- out$scale
	out$df.prior <- out$df2
	out$df2 <- out$scale <- NULL
	df.total <- df.residual + out$df.prior
	if(is.null(out$df.prior) || is.na(out$df.prior)) stop("Could not estimate prior df")
	if(out$df.prior == Inf)
		out$s2.post <- rep(out$s2.prior,length(sigma))
	else
		out$s2.post <- (ifelse(df.residual==0, 0, df.residual*sigma^2) + out$df.prior*out$s2.prior) / df.total
	out$t <- coefficients / stdev.unscaled / sqrt(out$s2.post)
	out$p.value <- 2*pt(-abs(out$t),df=df.total)

#	B-statistic
	if(is.null(std.coef)) {
		varpriorlim <- 10/out$s2.prior
		out$var.prior <- tmixture.matrix(out$t,stdev.unscaled,df.total,proportion,varpriorlim)
		out$var.prior[ is.na(out$var.prior) ] <- 1/out$s2.prior
		out$var.prior <- pmax(out$var.prior, 0.1/out$s2.prior)
	} else
		out$var.prior <- rep(std.coef^2/out$s2.prior, NCOL(out$t))
	r <- rep(1,NROW(out$t)) %o% out$var.prior
	r <- (stdev.unscaled^2+r) / stdev.unscaled^2
	t2 <- out$t^2
	if(out$df.prior > 10^6)
		kernel <- t2*(1-1/r)/2
	else
		kernel <- (1+df.total)/2*log((t2+df.total) / (t2/r+df.total))
	out$lods <- drop( log(proportion/(1-proportion))-log(r)/2+kernel )
	out
}

tmixture.matrix <- function(tstat,stdev.unscaled,df,proportion,c0lim=NULL) {
#	Estimate the prior variance of the coefficients for DE genes
#	Gordon Smyth
#	18 Nov 2002

	tstat <- as.matrix(tstat)
	stdev.unscaled <- as.matrix(stdev.unscaled)
	if(any(dim(tstat) != dim(stdev.unscaled))) stop("Dims of tstat and stdev.unscaled don't match")
	ncoef <- ncol(tstat)
	c0 <- rep(0,ncoef)
	for (j in 1:ncoef) c0[j] <- tmixture.vector(tstat[,j],stdev.unscaled[,j],df,proportion,c0lim)	
	c0
}

tmixture.vector <- function(tstat,stdev.unscaled,df,proportion,c0lim=NULL) {
#	Estimate scale factor in mixture of two t-distributions
#	tstat is assumed to follow (c0+c1)/c1*t(df) with probability proportion and t(df) otherwise
#	Gordon Smyth
#	18 Nov 2002

	if(any(is.na(tstat))) {
		sel <- !is.na(tstat)
		tstat <- tstat[sel]
		stdev.unscaled <- stdev.unscaled[sel]
		df <- df[sel]
	}
	ngenes <- length(tstat)
	ntarget <- ceiling(proportion/2*ngenes)
	if(ntarget < 1) return(NA)
	tstat <- abs(tstat)
	ttarget <- quantile(tstat,(ngenes-ntarget)/(ngenes-1))
	top <- (tstat >= ttarget)
	tstat <- tstat[top]
	c1 <- stdev.unscaled[top]^2
	df <- df[top]
	r <- ntarget-rank(tstat)+1
	p0 <- pt(-tstat,df=df)
	ptarget <- ( (r-0.5)/2/ngenes - (1-proportion)*p0 ) / proportion
	pos <- ptarget > p0
	c0 <- rep(0,ntarget)
	if(any(pos)) {
		qtarget <- qt(ptarget[pos],df=df[pos])
		c0[pos] <- c1[pos]*((tstat[pos]/qtarget)^2-1)
	}
	if(!is.null(c0lim)) c0 <- pmin(c0,c0lim)
	mean(c0)
}

fitFDist <- function(x,df1) {
#	Moment estimation of the parameters of a scaled F-distribution
#	The first degrees of freedom is given
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 Apr 2003.

#	Remove missing or infinite values and zero degrees of freedom
	o <- is.finite(x) & is.finite(df1) & (x > 0) & (df1 > 0)
	if(any(!o)) {
		x <- x[o]
		df1 <- df1[o]
	}
	n <- length(x)
	if(n==0) return(list(scale=NA,df2=NA))
	
#	Better to work on with log(F)
	z <- log(x)
	e <- z-digamma(df1/2)+log(df1/2)
	emean <- mean(e)
	evar <- mean(n/(n-1)*(e-emean)^2-trigamma(df1/2))
	if(evar > 0) {
		df2 <- 2*trigammaInverse(evar)
		s20 <- exp(emean+digamma(df2/2)-log(df2/2))
	} else {
		df2 <- Inf
		s20 <- exp(emean)
	}
	list(scale=s20,df2=df2)
}

trigammaInverse <- function(x) {
#	Solve trigamma(y) = x for y
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 April 2003.

#	Non-numeric or zero length input
	if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
	if(length(x)==0) return(numeric(0))

#	Treat out-of-range values as special cases
	omit <- is.na(x)
	if(any(omit)) {
		y <- x
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 0)
	if(any(omit)) {
		y <- x
		y[omit] <- NaN
		warning("NaNs produced")
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x > 1e7)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/sqrt(x[omit])
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 1e-6)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/x[omit]
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}

#	Newton's method
#	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
#	so iteration to solve 1/x = 1/trigamma is monotonically convergent
	y <- 0.5+1/x
	iter <- 0
	repeat {
		iter <- iter+1
		tri <- trigamma(y)
		dif <- tri*(1-tri/x)/tetragamma(y)
		y <- y+dif
		if(max(-dif/y) < 1e-8) break
		if(iter > 50) {
			warning("Iteration limit exceeded")
			break
		}
	}
	y
}

qqt <- function(y,df=Inf,ylim=range(y),main="Student's t Q-Q Plot",xlab="Theoretical Quantiles",ylab="Sample Quantiles",plot.it=TRUE,...)
{
#	Student's t probability plot
#	Gordon Smyth
#	3 Oct 2002

    y <- y[!is.na(y)]
    if(0 == (n <- length(y))) stop("y is empty")
    x <- qt(ppoints(n),df=df)[order(order(y))]
    if (plot.it) plot(x,y,main=main,xlab=xlab,ylab=ylab,ylim=ylim,...)
    invisible(list(x=x,y=y))
}

topTable <- function(fit,coef=1,number=10,genelist=NULL,adjust.method="holm",sort.by="B") {
#	Summary table of top genes, object-orientated version
#	Gordon Smyth
#	4 August 2003

	if(!missing(genelist)) fit$genes <- genelist
	toptable(fit=fit[c("coefficients","stdev.unscaled")],
		coef=coef,
		number=number,
		genelist=fit$genes,
		A=fit$Amean,
		eb=fit[c("t","p.value","lods")],
		adjust.method=adjust.method,
		sort.by=sort.by)
}

toptable <- function(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="holm",sort.by="B",...) {
#	Summary table of top genes
#	Gordon Smyth
#	21 Nov 2002. Last revised 18 Sep 2003.

	if(is.null(eb)) {
		fit$coefficients <- as.matrix(fit$coef)[,coef]
		fit$stdev.unscaled <- as.matrix(fit$stdev)[,coef]
		eb <- ebayes(fit,...)
		coef <- 1
	}
	M <- as.matrix(fit$coef)[,coef]
	if(is.null(A)) {
		if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
	} else {
		if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
	}
	tstat <- as.matrix(eb$t)[,coef]
	P.Value <- as.matrix(eb$p)[,coef]
	B <- as.matrix(eb$lods)[,coef]
	ord <- switch(sort.by,
		M=order(abs(M),decreasing=TRUE),
		A=order(A,decreasing=TRUE),
		P=order(P.Value,decreasing=FALSE),
		p=order(P.Value,decreasing=FALSE),
		T=order(abs(tstat),decreasing=TRUE),
		t=order(abs(tstat),decreasing=TRUE),
		B=order(B,decreasing=TRUE),order(B,decreasing=TRUE))
	top <- ord[1:number]
	P.Value <- p.adjust(P.Value,method=adjust.method)
	if(is.null(genelist))
		tab <- data.frame(M=M[top])
	else if(is.null(dim(genelist)))
		tab <- data.frame(Name=I(genelist[top]),M=M[top])
	else
		tab <- data.frame(genelist[top,],M=M[top])
	if(!is.null(A)) tab <- data.frame(tab,A=A[top])
	tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],B=B[top])
	rownames(tab) <- as.character(1:length(M))[top]
	tab
}

