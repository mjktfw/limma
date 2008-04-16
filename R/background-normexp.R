#  BACKGROUND-NORMEXP.R

#  NORMAL + EXPONENTIAL ADAPTIVE MODEL

normexp.signal <- function(par,x)
#	Expected value of signal given foreground in normal + exponential model
#	Gordon Smyth
#	24 Aug 2002. Last modified 26 December 2005.
{
	mu <- par[1]
	sigma <- exp(par[2])
	alpha <- exp(par[3])
#	cat(c(mu,sigma,alpha),"\n")
	if(alpha <= 0) stop("alpha must be positive")
	if(sigma <= 0) stop("sigma must be positive")
	mu.sf <- x-mu-sigma^2/alpha
	signal <- mu.sf + sigma^2 * exp(dnorm(0,mean=mu.sf,sd=sigma,log=TRUE) - pnorm(0,mean=mu.sf,sd=sigma,lower.tail=FALSE,log=TRUE))
	o <- !is.na(signal)
	if(any(signal[o]<0)) {
		warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensity to small value")
		signal[o] <- pmax(signal[o],1e-6)
	}
	signal
}

normexp.fit <- function(x, method="saddle", n.pts=NULL, trace=FALSE)
#	Estimate parameters of normal+exponential convolution model
#	Gordon Smyth and Jeremy Silver
#	24 Aug 2002. Last modified 22 July 2007.
{
	isna <- is.na(x)
	if(any(isna)) x <- x[!isna]
	if(length(x)<4) stop("Not enough data: need at least 4 non-missing corrected intensities")

	method <- match.arg(method,c("saddle","neldermead","bfgs","rma","mcgee"))

	if(method=="rma") {
		require(affy)
		Results <- bg.parameters(x)
		return(list(par=c(Results$mu,log(Results$sigma),-log(Results$alpha))))
	}

	if(method=="mcgee") {
		Results <- .bg.parameters.rma75(x)
		return(list(par=c(Results$mu,log(Results$sigma),-log(Results$alpha))))
	}

#	Starting values for parameters mu, alpha and sigma
	q <- quantile(x, c(0,0.05,0.1,1), na.rm = TRUE, names = FALSE)
	if(q[1]==q[4]) return(list(par=c(q[1],0,0),m2loglik=NA,convergence=0))
	if(q[2] > q[1]) {
		mu <- q[2]
	} else {
		if(q[3] > q[1]) {
			mu <- q[3]
		} else {
			mu <- q[1] + 0.05*(q[4]-q[1])
		}
	}
	sigma <- sqrt(mean((x[x<mu]-mu)^2, na.rm = TRUE))
	alpha <- mean(x,na.rm = TRUE) - mu
	if(alpha <= 0) alpha <- 1e-6
#	if(trace) cat("Starting values\n",mu,sigma,alpha,"\n")

#	Use a maximum of n.pts points for the fit
	if(!is.null(n.pts)) if(n.pts >= 4 & n.pts < length(x)) {
		a <- 0.5
		x <- quantile(x,((1:n.pts)-a)/n.pts,type=5)
	}

	Results <- switch(method,
		"saddle"= optim(par=c(mu,log(sigma),log(alpha)), fn=normexp.m2loglik.saddle, control=list(trace=as.integer(trace)), x=x),
		"neldermean"= optim(par=c(mu,log(sigma),log(alpha)), fn=normexp.m2loglik, control=list(trace=as.integer(trace)), x=x),
		"bfgs"= optim(par=c(mu,log(sigma),log(alpha)), fn=normexp.m2loglik, gr=normexp.grad, method=c("BFGS"), control=list(trace=as.integer(trace)), x=x)
	)
	if(trace) cat("normexp par:",Results$par[1],exp(Results$par[2]),exp(Results$par[3]),"\n")
	list(par=Results$par,m2loglik=Results$value,convergence=Results$convergence)
}

normexp.grad <- function(par,x)
#	Gradient of norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver.
#	21 Jan 2005. Last modified 26 December 2005.
{
	mu <- par[1]
	logsigma <- par[2]
	logalpha <- par[3]
	e <- x-mu
	mu.sf <- e - exp(2*logsigma - logalpha)

	dlogdbeta <- -2 * sum(exp(-logalpha) - exp(dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))
	dlogdlogsigma <- -2 * sum(exp(2*logsigma - 2*logalpha) - (2*exp(2*logsigma - logalpha) + mu.sf)*exp( dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))
	dlogdlogalpha <- -2 * sum(-1 + e * exp( -logalpha) - exp(2*logsigma - 2*logalpha) + exp(2*logsigma - logalpha+dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))

	c(dlogdbeta,dlogdlogsigma,dlogdlogalpha)
}

dnormexp <- function(x,normmean,normsd,expmean,log=FALSE)
#	Density of normal+exponential convolution
#	Gordon Smyth and Jeremy Silver
#	24 Aug 2002. Last modified 26 Dec 2005.
{
	e <- (x-normmean)/normsd
	a2 <- expmean/normsd
	mu.sf <- e-1/a2
	logf <- -log(expmean) + (0.5/a2-e)/a2 + pnorm(0,mu.sf,1,lower.tail=FALSE,log.p=TRUE)
	if(log) logf else exp(logf) 
}

dnormexp.saddle <- function(x,normmean,normsd,expmean,log=FALSE,secondorder=TRUE)
#	Saddlepoint approximation to normal+exponential density
#	Gordon Smyth
#	26 December 2005.  Revised 19 July 2007.

#	The approx is nearly exact at the low end, and is under about 0.002 at the high end:
#	> dexp(64000,rate=1/15000,log=TRUE)-dnormexp.saddle(64000,0,1e-8,15000,log=TRUE)
#	[1] 0.002271867
{
	mu <- normmean
	sigma <- normsd
	alpha <- expmean
	theta <- .findSaddleTheta(x,mu,sigma,alpha)
	k <- mu*theta+0.5*sigma^2*theta^2-log(1-alpha*theta)
	k2 <- sigma^2+alpha^2/(1-alpha*theta)^2
	logf <- -0.5*log(2*pi*k2)-x*theta+k
	if(secondorder) {
		k3 <- 2*alpha^3/(1-alpha*theta)^3
		k4 <- 6*alpha^4/(1-alpha*theta)^4
		logf <- logf+1/8*k4/k2^2-5/24*k3^2/k2^3
	}
	if(log) logf else exp(logf)
}

.findSaddleTheta <- function(x,normmean,normsd,expmean,trace=FALSE)
#	Find the canonical parameter (theta) for the saddle-point approximation
#	Gordon Smyth
#	22 July 2007
{
	mu <- normmean
	sigma <- normsd
	alpha <- expmean
	e <- x-mu

#	Sigma small approximation
	upperbound1 <- pmax(0,(e-alpha)/alpha/abs(e))

#	alpha small approximation
	upperbound2 <- e/sigma^2
	upperbound <- pmin(upperbound1,upperbound2)

#	Solve quadratic approximation
#	Theoretically exact, but subject to subtractive cancellation
	c2 <- sigma^2*alpha
	c1 <- -sigma^2-e*alpha
	c0 <- -alpha+e
	theta.quadratic <- (-c1-sqrt(c1^2-4*c0*c2))/2/c2

#	Globally convergence Newton iteration
	theta <- pmin(theta.quadratic,upperbound)
	i <- 0
	repeat {
		i <- i+1
		dK <- mu+sigma^2*theta+alpha/(1-alpha*theta)
		if(trace) cat("max dev",max(abs(dK-x)),"\n")
		ddK <- sigma^2+(alpha/(1-alpha*theta))^2
		delta <- (x-dK)/ddK
		theta <- theta + delta
		if(i==1) theta <- pmin(theta,upperbound)
		if(all(abs(delta) < 1e-10)) break
		if(i > 50) break
	}
	if(trace) {
		dK <- mu+sigma^2*theta+alpha/(1-alpha*theta)
		plot(x,dK-x)
	}
	theta
}

normexp.m2loglik.saddle <- function(par,x)
#	Minus twice the norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver and Gordon Smyth
#	24 Aug 2002. Last modified 26 December 2005.
{
	-2*sum(dnormexp.saddle(x,normmean=par[1],normsd=exp(par[2]),expmean=exp(par[3]),log=TRUE,second=TRUE))
}

normexp.m2loglik <- function(par,x)
#	Minus twice the norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver and Gordon Smyth
#	24 Aug 2002. Last modified 26 December 2005.
{
	-2*sum(dnormexp(x,normmean=par[1],normsd=exp(par[2]),expmean=exp(par[3]),log=TRUE))
}

.bg.parameters.rma75 <- function(pm,n.pts = 2^14)
#	Estimate estimate normexp parameters
#	This is extracted from the RMA-75 function of
#	McGee, M. and Chen, Z. (2006). Parameter estimation for the
#	exponential-normal convolution model for background correction
#	of Affymetrix GeneChip data.
#	Stat Appl Genet Mol Biol, 5(1), Article 24.
{
##	mu-correction function
	mu.est.correct <- function(m,s,a) { 
		f <- function(x) (dnorm(x-s*a)-s*a*(pnorm(x-s*a)+pnorm(m/s+s*a)-1))
		t <- uniroot(f, c(-5, 10), tol = 1e-12)$root
		t <- m-s*t
		return(t)
	}

##	getting mode function
	max.density <- function(x, n.pts) {
		aux <- density(x, kernel = "epanechnikov", n = n.pts, na.rm = TRUE)
		aux$x[order(-aux$y)[1]]
	}

	pmbg <- max.density(pm, n.pts)
	bg.data <- pm[pm < pmbg]		   
	pmbg <- max.density(bg.data, n.pts) 
	mubg <- pmbg   ## the mode
	bg.data <- pm[pm < pmbg]
	bg.data <- bg.data - pmbg
	bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2) ## estimate sigma
	sig.data<-pm[pm > pmbg]
	sig.data <- sig.data - pmbg   
	q75 <- 0.75
	alpha3 <- -(quantile(pm,q75)-pmbg)/log(1-q75) ## 75th quantile estimation

##	mode-correction
	mu3 <- mu.est.correct(m=mubg,s=bgsd,a=1/alpha3)
	mu3 <- (mu3+mubg)/2  ## take ave
	bg.data3<- pm[pm < mu3] 
	bg.data3 <- bg.data3 - mu3
	bgsd3 <- sqrt(sum(bg.data3^2)/(length(bg.data3) - 1)) * sqrt(2)
	sig.data3 <- pm[pm > mu3]
	alpha3<- -(quantile(pm,q75)-mu3)/log(1-q75)
	list(alpha = 1/alpha3, mu = mu3, sigma = bgsd3)
}
