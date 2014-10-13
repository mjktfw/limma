predFCm <- function(fit,coef=2,prob=TRUE,VarRel=NULL)
# Belinda Phipson 29 May 2012. Updated 8 January 2013.
{
 p <- 1-propTrueNull(fit$p.value[,coef],method="lfdr")
 if(p==0) p<-1e-8
 if(length(fit$s2.prior)==1) trend<-FALSE else trend<-TRUE
 fit <- eBayes(fit,proportion = p,trend=trend)
 v <- fit$cov.coefficients[coef,coef]

 if(is.null(VarRel)) VarRel <- "Independent"

 if(VarRel=="Independent"){
	v0 <- fitGammaIntercept(fit$coeff[,coef]^2,offset=v*fit$s2.post)
	if(v0<0) v0<-1e-8
	pfc <- fit$coeff[,coef]*v0/(v0+v*fit$s2.post)
	if(prob){
		A <- p/(1-p)
		B <- (v*fit$s2.post/(v*fit$s2.post+v0))^0.5
		C <- exp(fit$coeff[,coef]^2*v0/(2*v^2*fit$s2.post^2+2*v*v0*fit$s2.post))
		lods <- log(A*B*C)
		probDE <- exp(lods)/(1+exp(lods))
		probDE[lods>700] <- 1
		pfc <- pfc*probDE
		}
	}
 else if(VarRel=="Increasing"){
	b2 <- fit$coeff[,coef]^2/fit$s2.post
        v0 <- fitGammaIntercept(b2,offset=v)
	if(v0<0) v0<-1e-8
	pfc <- fit$coeff[,coef]*v0/(v0+v)
	if(prob){
		A <- p/(1-p)
		B <- (v/(v+v0))^0.5
		C <- exp(fit$coeff[,coef]^2*v0/(2*v^2*fit$s2.post+2*v*v0*fit$s2.post))
		lods <- log(A*B*C)
		probDE <- exp(lods)/(1+exp(lods))
		probDE[lods>700] <- 1
		pfc <- pfc*probDE
		}
	}
 else stop("Invalid VarRel, please select either Independent or Increasing")
 pfc
}
