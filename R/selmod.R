##  selmod.R

selectModel <- function(y, designlist, criterion="bic", var.prior=NULL, df.prior=0, ...)
# y is a data matrix to be fitted, with rows as genes and columns as arrays
# designList is a list of design matrices to be compared
# the function returns
# BIC values for each design and
# the prefered model calculated as the minimum BIC from
# the compared models for each gene separately
# written 17/7/08 Alicia Oshlack
# Last revised 28 Aug 2008, Gordon Smyth
{
	y <- as.matrix(y)
	narrays <- ncol(y)
	nmodels <- length(designlist)
	models <- names(designlist)
	if(is.null(models)) models <- as.character(1:nmodels)
	if(df.prior>0 & is.null(var.prior)) stop("var.prior must be set")
	if(df.prior==0) var.prior <- 0
	criterion <- match.arg(criterion,c("bic","aic","lik"))
	ntotal <- df.prior+narrays
	penalty <- switch(criterion,bic=log(ntotal),aic=2,lik=0)
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

