#  removeBatchEffect.R

#  A refinement would be to empirical Bayes shrink
#  the batch effects before subtracting them.

removeBatchEffect <- function(x,batch,batch2=NULL,design=matrix(1,ncol(x),1))
#  Remove batch effects from matrix of expression data
#  Carolyn de Graaf and Gordon Smyth
#  Created 1 Aug 2008. Last revised 31 Aug 2009.
{
	x <- as.matrix(x)
	batch <- as.factor(batch)
	contrasts(batch) <- contr.sum(levels(batch))
	if(is.null(batch2)) {
		X <- model.matrix(~batch)[,-1,drop=FALSE]
	} else {
		batch2 <- as.factor(batch2)
		contrasts(batch2) <- contr.sum(levels(batch2))
		X <- model.matrix(~batch+batch2)[,-1,drop=FALSE]
	}
	X <- qr.resid(qr(design),X)
	t(qr.resid(qr(X),t(x)))
}

