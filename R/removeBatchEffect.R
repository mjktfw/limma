#  removeBatchEffect.R

removeBatchEffect <- function(x,batch,design=NULL)
#  Remove batch effects from matrix of expression data
#  Carolyn de Graaf and Gordon Smyth
#  1 Aug 2008. Last revised 15 Oct 2008.
{
	x <- as.matrix(x)
	batch <- as.factor(batch)
	contrasts(batch) <- contr.sum(levels(batch))
	X <- model.matrix(~batch)[,-1,drop=FALSE]
	if(!is.null(design)) X <- qr.resid(qr(design),X)
	qrX <- qr(X)
	t(qr.resid(qrX,t(x)))
}

