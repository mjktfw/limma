designMatrix <- function(targets, ref) {
#	Design matrix for two-color experiments
#	'targets' is matrix or data.frame with columns Cy3 and Cy5
#	Gordon Smyth
#	25 June 2003

	tar <- targets[,c("Cy3","Cy5")]
	tar <- as.vector(t(as.matrix(tar)))
	lev <- unique(tar)
	treatments <- setdiff(lev,ref)
	lev <- c(ref,treatments)
	tar <- factor(tar,levels=lev)
	n <- length(tar)
	col <- factor(rep(c(1,2),length=n))
	contrasts(col) <- matrix(c(-1,1),2,1)
	X <- model.matrix(~-1+tar*col)
	keeprows <- X[,1]==0
	keepcols <- attr(X,"assign")==3
	design <- X[keeprows,keepcols]
	rownames(design) <- targets$SlideNumber
	colnames(design) <- treatments
	design
}
