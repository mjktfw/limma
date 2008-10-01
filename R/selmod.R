##  selmod.R

selectModel <- function(y, designlist, criterion="aic", var.prior=NULL, df.prior=0, ...)
# y is a data matrix to be fitted, with rows as genes and columns as arrays.
# designlist is a list of design matrices to be compared.
# The function returns AIC or BIC values for each design
# and the prefered model for each gene with minimim criterion value.
# Written 17/7/08 Alicia Oshlack
# Last revised 1 Oct 2008, Gordon Smyth
{
	y <- as.matrix(y)
	narrays <- ncol(y)
	nmodels <- length(designlist)
	models <- names(designlist)
	if(is.null(models)) models <- as.character(1:nmodels)
	if(df.prior>0 & is.null(var.prior)) stop("var.prior must be set")
	if(df.prior==0) var.prior <- 0
	criterion <- match.arg(criterion,c("aic","bic","lik"))
	ntotal <- df.prior+narrays
	penalty <- switch(criterion,bic=log(narrays),aic=2,lik=0)
	for(i in 1:nmodels) {
		fit <- lmFit(y, designlist[[i]], ...)
		npar <- narrays-fit$df.residual
		s2.post <- (df.prior*var.prior+fit$df.residual*fit$sigma^2)/ntotal
		if(i==1) IC <- matrix(nrow=nrow(fit),ncol=nmodels,dimnames=list(Probes=rownames(fit),Models=models))
		IC[,i] <- ntotal*log(s2.post)+npar*penalty
	}
	pref <- factor(apply(IC,1,which.min),levels=1:nmodels,labels=models)
	list(IC=IC,pref=pref,criterion=criterion)
}

