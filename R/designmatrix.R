designMatrix <- function(targets, ref) {
#	Design matrix for two-color experiments
#	'targets' is matrix or data.frame with columns Cy3 and Cy5
#	Gordon Smyth
#	25 June 2003.  Last modified 29 June 2003.

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
	design <- X[keeprows,keepcols,drop=FALSE]
	rownames(design) <- rownames(targets)
	colnames(design) <- treatments
	design
}

makeContrasts <- function(..., levels) {
#	Construct matrix of custom contrasts
#	Gordon Smyth
#	30 June 2003.

	if(is.factor(levels)) levels <- levels(levels)
	if(is.matrix(levels)) levels <- colnames(levels)
	levels <- make.names(levels)
	n <- length(levels)
	if(n < 1) stop("No levels to construct contrasts from")
	indicator <- function(i,n) {
		out <- rep(0,n)
		out[i] <- 1
		out
	}
	for (i in 1:n) assign(levels[i], indicator(i,n))
	e <- substitute(list(...))
	ne <- length(e)
	cm <- matrix(0,n,ne-1)
	rownames(cm) <- levels
	if(ne < 2) return(cm)
	colnames(cm) <- as.character(e)[2:ne]
	for (j in 1:(ne-1)) {
		ej <- e[[j+1]]
		if(is.character(ej)) ej <- parse(text=ej)
		cm[,j] <- eval(ej)
	}
	cm
}
