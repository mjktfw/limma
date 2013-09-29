#  SUBSET DATA SETS

assign("[.RGList",
function(object, i, j, ...) {
#  Subsetting for RGList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 22 December 2005.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	oc <- names(object$other)
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
			for(k in oc) object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$R <- object$R[i,,drop=FALSE]
			object$G <- object$G[i,,drop=FALSE]
			object$Rb <- object$Rb[i,,drop=FALSE]
			object$Gb <- object$Gb[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			for(k in oc) object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
		} else {
			object$R <- object$R[i,j,drop=FALSE]
			object$G <- object$G[i,j,drop=FALSE]
			object$Rb <- object$Rb[i,j,drop=FALSE]
			object$Gb <- object$Gb[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			for(k in oc) object$other[[k]] <- object$other[[k]][i,j,drop=FALSE]
		}
	}
	object
})

assign("[.MAList",
function(object, i, j, ...) {
#  Subsetting for MAList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 22 Dec 2005.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	other <- names(object$other)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$M <- object$M[,j,drop=FALSE]
			object$A <- object$A[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- as.matrix(object$design)[j,,drop=FALSE]
				if(!is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)
			}
			for(a in other) object$other[[a]] <- object$other[[a]][,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$M <- object$M[i,,drop=FALSE]
			object$A <- object$A[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			for(a in other) object$other[[a]] <- object$other[[a]][i,,drop=FALSE]
		} else {
			object$M <- object$M[i,j,drop=FALSE]
			object$A <- object$A[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- as.matrix(object$design)[j,,drop=FALSE]
				if(!is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)
			}
			for(a in other) object$other[[a]] <- object$other[[a]][i,j,drop=FALSE]
		}
	}
	object
})

assign("[.EList",
function(object, i, j, ...) {
#  Subsetting for EList objects
#  Gordon Smyth
#  23 February 2009.  Last modified 21 October 2010.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	other <- names(object$other)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$E <- object$E[,j,drop=FALSE]
			object$Eb <- object$Eb[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- as.matrix(object$design)[j,,drop=FALSE]
				if(!is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)
			}
			for(a in other) object$other[[a]] <- object$other[[a]][,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$E <- object$E[i,,drop=FALSE]
			object$Eb <- object$Eb[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			for(a in other) object$other[[a]] <- object$other[[a]][i,,drop=FALSE]
		} else {
			object$E <- object$E[i,j,drop=FALSE]
			object$Eb <- object$Eb[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- as.matrix(object$design)[j,,drop=FALSE]
				if(!is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)
			}
			for(a in other) object$other[[a]] <- object$other[[a]][i,j,drop=FALSE]
		}
	}
	object
})

assign("[.EListRaw", get("[.EList"))

assign("[.MArrayLM",
function(object, i, j, ...)
#  Subsetting for MArrayLM objects
#  Gordon Smyth
#  26 April 2005. Last modified 28 September 2013.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(!is.null(object$coefficients)) object$coefficients <- as.matrix(object$coefficients)
	if(!is.null(object$stdev.unscaled)) object$stdev.unscaled <- as.matrix(object$stdev.unscaled)
	if(!is.null(object$weights)) object$weights <- as.matrix(object$weights)
	if(!is.null(object$p.value)) object$p.value <- as.matrix(object$p.value)
	if(!is.null(object$lods)) object$lods <- as.matrix(object$lods)
	if(!is.null(object$targets)) object$targets <- as.data.frame(object$targets)
	if(!is.null(object$cov.coefficients)) object$cov.coefficients <- as.matrix(object$cov.coefficients)
	if(!is.null(object$contrasts)) object$contrasts <- as.matrix(object$contrasts)
	if(is.null(object$contrasts) && !is.null(object$coefficients)) {
		object$contrasts <- diag(ncol(object$coefficients))
		rownames(object$contrasts) <- colnames(object$contrasts) <- colnames(object$coefficients)
	}
	if(!is.null(object$genes)) object$genes <- as.data.frame(object$genes)
	if(missing(i)) {
		if(missing(j))
			return(object)
		else {
			object$coefficients <- object$coefficients[,j,drop=FALSE]
			object$stdev.unscaled <- object$stdev.unscaled[,j,drop=FALSE]
			object$t <- object$t[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$p.value <- object$p.value[,j,drop=FALSE]
			object$lods <- object$lods[,j,drop=FALSE]
			object$cov.coefficients <- object$cov.coefficients[j,j,drop=FALSE]
			object$contrasts <- object$contrasts[,j,drop=FALSE]
			object$var.prior <- object$var.prior[j]
		}
	} else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$coefficients <- object$coefficients[i,,drop=FALSE]
			object$stdev.unscaled <- object$stdev.unscaled[i,,drop=FALSE]
			object$t <- object$t[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$p.value <- object$p.value[i,,drop=FALSE]
			object$lods <- object$lods[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
		} else {
			object$coefficients <- object$coefficients[i,j,drop=FALSE]
			object$stdev.unscaled <- object$stdev.unscaled[i,j,drop=FALSE]
			object$t <- object$t[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$p.value <- object$p.value[i,j,drop=FALSE]
			object$lods <- object$lods[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$cov.coefficients <- object$cov.coefficients[j,j,drop=FALSE]
			object$contrasts <- object$contrasts[,j,drop=FALSE]
			object$var.prior <- object$var.prior[j]
		}
		object$df.residual <- object$df.residual[i]
		if(length(object$df.prior)>1) object$df.prior <- object$df.prior[i]
		object$df.total <- object$df.total[i]
		object$sigma <- object$sigma[i]
		object$s2.post <- object$s2.post[i]
		object$Amean <- object$Amean[i]
	}
	if(!is.null(object$F))
		if(missing(j)) {
			object$F <- object$F[i]
			object$F.p.value <- object$F.p.value[i]
		} else {
			F.stat <- classifyTestsF(object,fstat.only=TRUE)
			object$F <- as.vector(F.stat)
			df1 <- attr(F.stat,"df1")
			df2 <- attr(F.stat,"df2")
			if (df2[1] > 1e+06) 
				object$F.p.value <- pchisq(df1*object$F,df1,lower.tail=FALSE)
			else
				object$F.p.value <- pf(object$F,df1,df2,lower.tail=FALSE)
		}
	object
})

