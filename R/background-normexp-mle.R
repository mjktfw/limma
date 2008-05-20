normexp.fit.C <- function(x, method="saddle", n.pts=NULL, trace=FALSE)
#	Estimate parameters of normal+exponential convolution model
#	Pure R Version Gordon Smyth 24 Aug 2002.
#	Version with C by Jeremy Silver 29 Oct 2007.
#	Last modified 9 May 2008.
{
	isna <- is.na(x)
	if(any(isna)) x <- x[!isna]
	if(length(x)<4) stop("Not enough data: need at least 4 non-missing corrected intensities")

	method <- match.arg(method,c("saddle","neldermead","bfgs","rma","mcgee","nlminb","nlminblog"))

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
	s2 <- sigma^2;
	alpha <- mean(x,na.rm = TRUE) - mu
	if(alpha <= 0) alpha <- 1e-6
#	if(trace) cat("Starting values\n",mu,sigma,alpha,"\n")

#	Use a maximum of n.pts points for the fit
	if(!is.null(n.pts)) if(n.pts >= 4 & n.pts < length(x)) {
		a <- 0.5
		x <- quantile(x,((1:n.pts)-a)/n.pts,type=5)
	}

	Results <- switch(method,
		"saddle"= .nelderMeanSaddleInC(Pars = c(mu, log(sigma), log(alpha)),X = x),
		"neldermead"= optim(par=c(mu,log(sigma),log(alpha)), fn=normexp.m2loglik, control=list(trace=as.integer(trace)), x=x),
		"bfgs"= optim(par=c(mu,log(sigma),log(alpha)), fn=normexp.m2loglik, gr=normexp.grad, method=c("BFGS"), control=list(trace=as.integer(trace)), x=x),
		"nlminb" = nlminb(start = c(mu,alpha,s2),objective = .normexp.m2sumloglik.alternate,gradient = .normexp.gm2sumloglik.alternate, hessian = .normexp.hm2sumloglik.alternate, f = x,bg = rep(0,length(x)),control = list(trace = as.integer(trace)),scale = median(c(mu,alpha,s2))/c(mu,alpha,s2)),
		"nlminblog" = .nlminbLog(Pars = c(mu, log(sigma), log(alpha)),X = x)
	)
	if(trace) cat("normexp par:",Results$par[1],exp(Results$par[2]),exp(Results$par[3]),"\n")
	if(method == "nlminb")
		return(list(par=c(Results$par[1],0.5*log(Results$par[2]),log(Results$par[3])),m2loglik=Results$objective,convergence=Results$convergence))
	else if(method == "nlminblog")
		return(list(par=c(Results$par[1],0.5*Results$par[2],Results$par[3]),m2loglik=Results$objective,convergence=Results$convergence))
	else
		return(list(par=Results$par,m2loglik=Results$value,convergence=Results$convergence))
}

.nelderMeanSaddleInC <- function(Pars,X)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  Out <- .C("fit_saddle_nelder_mead",
    parsIn = as.double(Pars), 
    X = as.double(X), 
    N = as.integer(length(X)), 
    fail = as.integer(0), 
    fncount = as.integer(0), 
    Fmin = as.double(0),
    PACKAGE="limma")
  list(par = Out$parsIn, value = Out$Fmin, convergence = Out$fail, fncount = Out$fncount)
}

.nlminbLog <- function(Pars,X)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  Out1 <- .C("fit_saddle_nelder_mead",
    parsIn = as.double(Pars), 
    X = as.double(X), 
    N = as.integer(length(X)), 
    fail = as.integer(0), 
    fncount = as.integer(0), 
    Fmin = as.double(0),
    PACKAGE="limma")[c(1,4:6)]
    
  Pars2 <- Out1$parsIn;
  Pars2[2] <- 2*Pars2[2];
  initLL <- .normexp.m2sumloglik.alternate3(Pars2,f = X,bg = rep(0,length(X)))

  Out2 <- nlminb(start = Pars2,objective = .normexp.m2sumloglik.alternate3,gradient = .normexp.gm2sumloglik.alternate3, hessian = .normexp.hm2sumloglik.alternate3, f = X,bg = rep(0,length(X)),scale = median(abs(Pars2))/abs(Pars2))
  if(Out2$objective < initLL){
    return(list(par = Out2$par, objective = Out2$objective, convergence = Out2$convergence))
  } else {
    return(list(par = Pars2, objective = initLL, convergence = Out1$fail))
  }
}  

.normexp.m2sumloglik.alternate <- function(theta,f,bg)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  bt <- theta[1]
  s2 <- theta[2]
  al <- theta[3]
  mu <- bt + bg
  e <- f - mu
  mu.sf <- e - s2/al

  -2*sum(-log(al) - e/al + 0.5*s2/(al^2) + pnorm(0,mu.sf,sqrt(s2),lower.tail = FALSE,log.p = TRUE))
}

.normexp.gm2sumloglik.alternate <- function(theta,f,bg)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  bt <- theta[1]
  s2 <- theta[2]
  al <- theta[3]
  mu <- bt + bg
  e <- f - mu
  mu.sf <- e - s2/al
  psionPsi <- dnorm(0,mu.sf,sqrt(s2))/pnorm(0,mu.sf,sqrt(s2),lower.tail = FALSE)

  dL.dbt <- sum(1/al - psionPsi)
  dL.ds2 <- sum(0.5/(al^2) - (1/al + 0.5*mu.sf/s2) * psionPsi)
  dL.dal <- sum(e/(al^2) - 1/al - s2/(al^3) + psionPsi*s2/(al^2))
  -2*c(dL.dbt,dL.ds2,dL.dal)

}

.normexp.hm2sumloglik.alternate <- function(theta,f,bg)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  bt <- theta[1]
  s2 <- theta[2]
  al <- theta[3]
  mu <- bt + bg
  e <- f - mu
  s2onal <- s2/al
  mu.sf <- e - s2onal
  psionPsi <- dnorm(0,mu.sf,sqrt(s2))/pnorm(0,mu.sf,sqrt(s2),lower.tail = FALSE)
  psionPsi2 <- psionPsi^2

  d2L.dbtdbt <- sum(-psionPsi2 - psionPsi*mu.sf/s2)
  d2L.dbtds2 <- sum( -0.5*(e + s2onal)*psionPsi2/s2 + 0.5*(-((e + s2onal)^2) + 2*s2onal*(e + s2onal) + s2)*psionPsi/(s2^2)) # OK TO 3 DP
  d2L.dbtdal <- sum( -al^-2 + s2onal*psionPsi2/al + mu.sf*psionPsi/(al^2))
  d2L.ds2ds2 <- sum( -(0.25/(s2^2))*((e + s2onal)^2)*psionPsi2 + psionPsi*(-e^3 + e*(3*al - e)*s2onal + (e + al)*(s2onal^2) + (s2onal^3))/(4*(s2^3)) )
  d2L.dalds2 <- sum( -1/(al^3) + (al^-2)*0.5*(psionPsi2*(e + s2onal) + (e^2 + s2 - (s2onal^2))*psionPsi/s2))
  d2L.daldal <- sum( (al^-2) - 2*e/(al^3) + 3*s2/(al^4) - psionPsi2*((s2^2)/(al^4)) - psionPsi*(mu.sf + 2*al)*((s2)/(al^4)))
  -2*rbind(
  c(d2L.dbtdbt,d2L.dbtds2,d2L.dbtdal),
  c(d2L.dbtds2,d2L.ds2ds2,d2L.dalds2),
  c(d2L.dbtdal,d2L.dalds2,d2L.daldal))
}

.normexp.m2sumloglik.alternate3 <- function(theta,f,bg)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  bt <- theta[1]
  s2 <- exp(theta[2])
  al <- exp(theta[3])
  mu <- bt + bg
  
  .C("normexp_m2sumloglik_alternate",
    mu = as.double(mu), 
    s2 = as.double(s2), 
    al = as.double(al), 
    n = as.integer(length(f)), 
    f = as.double(f), 
    m2LL = double(1),
    PACKAGE="limma"
  )$m2LL
}

.normexp.gm2sumloglik.alternate3 <- function(theta,f,bg)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  bt <- theta[1]
  s2 <- exp(theta[2])
  al <- exp(theta[3])
  mu <- bt + bg

  .C("normexp_gm2sumloglik_alternate",
    mu = as.double(mu), 
    s2 = as.double(s2), 
    al = as.double(al), 
    n = as.integer(length(f)), 
    f = as.double(f), 
    dm2LL = double(3),
    PACKAGE = "limma"
  )$dm2LL

}

.normexp.hm2sumloglik.alternate3 <- function(theta,f,bg)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 9 May 2008.
{
  bt <- theta[1]
  s2 <- exp(theta[2])
  al <- exp(theta[3])
  mu <- bt + bg

  matrix(.C("normexp_hm2sumloglik_alternate",
    mu = as.double(mu), 
    s2 = as.double(s2), 
    al = as.double(al), 
    n = as.integer(length(f)), 
    f = as.double(f), 
    d2m2LL = double(9),
    PACKAGE="limma"
  )$d2m2LL,3,3)
}
