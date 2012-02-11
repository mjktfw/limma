#  removeBatchEffect.R

#  A refinement would be to empirical Bayes shrink
#  the batch effects before subtracting them.

removeBatchEffect <- function(x,batch=NULL,batch2=NULL,covariates=NULL,design=matrix(1,ncol(x),1))
#  Remove batch effects from matrix of expression data
#  Gordon Smyth and Carolyn de Graaf
#  Created 1 Aug 2008. Last revised 4 Feb 2012.
{
	x <- as.matrix(x)
	if(is.null(batch) && is.null(batch2) && is.null(covariates)) return(x)
	if(!is.null(batch)) {
		batch <- as.factor(batch)
		contrasts(batch) <- contr.sum(levels(batch))
		batch <- model.matrix(~batch)[,-1,drop=FALSE]
	}
	if(!is.null(batch2)) {
		batch2 <- as.factor(batch2)
		contrasts(batch2) <- contr.sum(levels(batch2))
		batch2 <- model.matrix(~batch2)[,-1,drop=FALSE]
	}
	if(!is.null(covariates)) covariates <- as.matrix(covariates)
	X <- cbind(batch,batch2,covariates)
	X <- qr.resid(qr(design),X)
	t(qr.resid(qr(X),t(x)))
}

