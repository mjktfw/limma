#  removeBatchEffect.R

#  A refinement would be to empirical Bayes shrink
#  the batch effects before subtracting them.

removeBatchEffect <- function(x,batch,design=matrix(1,ncol(x),1))
#  Remove batch effects from matrix of expression data
#  Carolyn de Graaf and Gordon Smyth
#  1 Aug 2008. Last revised 31 Aug 2009.
{
	x <- as.matrix(x)
	batch <- as.factor(batch)
	contrasts(batch) <- contr.sum(levels(batch))
	X <- model.matrix(~batch)[,-1,drop=FALSE]
	X <- qr.resid(qr(design),X)
	qrX <- qr(X)
	t(qr.resid(qrX,t(x)))
}

