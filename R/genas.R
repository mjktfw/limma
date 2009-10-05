##  GENAS.R

.multTLogLikNull <- function(x,fit,m) 
#	Calculate the log-likelihood under the null hypothesis of no biological correlation
#	Belinda Phipson and Gordon Smyth
#	21 September 2009. Last modified 21 September 2009.
{
	d0<-fit$df.prior[1]
	d<-fit$df.residual[1]
	s<-fit$s2.post
	B<-fit$coefficients
	m<-m
	V<-fit$cov.coefficients
	a1<-x[1]
	a2<-x[2]
	V0<-matrix(c(exp(a1),0,0,exp(a2)),2,2)
	
	y1<-log(gamma(0.5*(m+d0+d)))
	y2<-(0.5*m)*log(pi*(d0+d))
	y3<-log(gamma(0.5*(d0+d)))

	First<-y1-y2-y3

	R<-chol(V0+V)
	R2<-chol(V)

	Second<-sum(log(diag(R2)))-sum(log(diag(R)))

	W<-backsolve(R,t(B),transpose=TRUE)
	Q<-colSums(W^2)

	Third<-0.5*(m+d0+d)*log(1+Q/(s*(d0+d)))

	sum(First+Second-Third)
}


.multTLogLik <- function(x,fit,m) 
#            Calculate the log-likelihood with biological correlation
#     	     Belinda Phipson and Gordon Smyth
#    	     21 September 2009. Last modified 21 September 2009.
{	
 	     d0<-fit$df.prior[1]
 	     d<-fit$df.residual[1]
 	     s<-fit$s2.post
 	     B<-fit$coefficients
 	     m<-m
 	     V<-fit$cov.coefficients
 	     a1<-x[1]
 	     a2<-x[2]
 	     b<-x[3]
 	     L<-matrix(c(1,b,0,1),2,2)
 	     D<-matrix(c(exp(a1),0,0,exp(a2)),2,2)
	
 	     V0<-L %*% D %*% t(L)
	
 	     y1<-log(gamma(0.5*(m+d0+d)))
 	     y2<-(0.5*m)*log(pi*(d0+d))
 	     y3<-log(gamma(0.5*(d0+d)))
	
 	     First<-y1-y2-y3
	
 	     R<-chol(V0+V)
 	     R2<-chol(V)

 	     Second<-sum(log(diag(R2)))-sum(log(diag(R)))

 	     W<-backsolve(R,t(B),transpose=TRUE)
 	     Q<-colSums(W^2)

 	     Third<-0.5*(m+d0+d)*log(1+Q/(s*(d0+d)))

 	     sum(First+Second-Third)
}

genas <- function(fit,coef=c(1,2))
#     Genuine association of gene expression profiles
#     Belinda Phipson and Gordon Smyth
#     21 September 2009. Last modified 24 September 2009.
{
      if(ncol(fit)>2) fit<-fit[,coef]
      x<-log(fit$var.prior)
      m<-2
	
      Q2<-optim(x,.multTLogLikNull,fit=fit,m=m,control=list(fnscale=-2))

      Q1<-optim(c(Q2$par[1],Q2$par[2],0),.multTLogLik,fit=fit,m=m,control=list(fnscale=-2))

      L<-matrix(c(1,Q1$par[3],0,1),2,2)
      D<-matrix(c(exp(Q1$par[1]),0,0,exp(Q1$par[2])),2,2)
      V0<-L%*%D%*%t(L)
      rho<-V0[2,1]/sqrt(V0[1,1]*V0[2,2])

      V<-fit$cov.coefficients
      rhotech<-V[2,1]/sqrt(V[1,1]*V[2,2])

      D<--2*(Q2$value-Q1$value)
      p.val<-pchisq(D,df=1,lower.tail=FALSE)

      list(technical.correlation=rhotech,covariance.matrix=V0,biological.correlation=rho,deviance=D,p.value=p.val)
}