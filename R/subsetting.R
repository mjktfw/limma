#  SUBSET AND COMBINE DATA SETS

assign("[.RGList",
function(object, i, j, ...) {
#  Subsetting for RGList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 5 July 2003.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$R <- object$R[,j,drop=FALSE]
			object$G <- object$G[,j,drop=FALSE]
			object$Rb <- object$Rb[,j,drop=FALSE]
			object$Gb <- object$Gb[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
		}
	else
		if(missing(j)) {
			object$R <- object$R[i,,drop=FALSE]
			object$G <- object$G[i,,drop=FALSE]
			object$Rb <- object$Rb[i,,drop=FALSE]
			object$Gb <- object$Gb[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
		} else {
			object$R <- object$R[i,j,drop=FALSE]
			object$G <- object$G[i,j,drop=FALSE]
			object$Rb <- object$Rb[i,j,drop=FALSE]
			object$Gb <- object$Gb[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
		}
	object
})

assign("[.MAList",
function(object, i, j, ...) {
#  Subsetting for MAList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 5 July 2003.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$M <- object$M[,j,drop=FALSE]
			object$A <- object$A[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- object$design[j,,drop=FALSE]
				warning("subsetting design matrix component of MAList may make it singular",call.=FALSE)
			}
		}
	else
		if(missing(j)) {
			object$M <- object$M[i,,drop=FALSE]
			object$A <- object$A[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
		} else {
			object$M <- object$M[i,j,drop=FALSE]
			object$A <- object$A[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- object$design[j,,drop=FALSE]
				warning("subsetting design matrix component of MAList may make it singular",call.=FALSE)
			}
		}
	object
})

cbind.RGList <- function(..., deparse.level=1) {
#  Combine MAList objects assuming same genelists
#  Gordon Smyth
#  27 June 2003

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- cbind(out$R,objects[[i]]$R)
		out$G <- cbind(out$G,objects[[i]]$G)
		out$Rb <- cbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- cbind(out$Gb,objects[[i]]$Gb)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
	}
	out
}

cbind.MAList <- function(..., deparse.level=1) {
#  Combine MAList objects assuming same genelists
#  Gordon Smyth
#  27 June 2003

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$M <- cbind(out$M,objects[[i]]$M)
		out$A <- cbind(out$A,objects[[i]]$A)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
	}
	out
}

makeUnique <- function(x) {
#  Add characters to the elements of a character vector to make all values unique
#  Gordon Smyth
#  10 April 2003

	x <- as.character(x)
	tab <- table(x)
	tab <- tab[tab>1]
	lentab <- length(tab)
	if(lentab > 0) {
		u <- names(tab)
		for (i in 1:lentab) {
			n <- tab[i]
			x[x==u[i]] <- paste(x[x==u[i]],formatC(1:n,width=1+floor(log(n,10)),flag="0"),sep="")
		}
	}
	x
}

if(!isGeneric("merge")) setGeneric("merge")

setMethod("merge", c("RGList","RGList"), definition=
function(x,y,z,...) {
#  Merge RGList y into x aligning by row names
#  Gordon Smyth
#  11 April 2003

	genes1 <- rownames(x$R)
	if(is.null(genes1)) genes1 <- rownames(x$G)
	genes2 <- rownames(y$R)
	if(is.null(genes2)) genes2 <- rownames(y$G)
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two RGLists have different elements")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	for (i in fields1) x[[i]] <- cbind(x[[i]],y[[i]][ord2,])
	x
})
