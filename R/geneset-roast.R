##  ROAST.R

setClass("Roast",
#  rotation gene set test
representation("list")
)

setMethod("show","Roast",
#  Di Wu, Gordon Smyth
#  14 May 2010.  Last modified 19 May 2010.
function(object) print(object$p.value)
)

roast <- function(iset=NULL,y,design,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999)
# Rotation gene set testing for linear models
# Gordon Smyth and Di Wu
# Created 24 Apr 2008. Revised 10 May 2012.
{
	if(is.null(iset)) iset <- rep(TRUE,nrow(y))
	y <- as.matrix(y)
	design <- as.matrix(design)
	if(!is.null(df.prior) && df.prior<0) stop("df.prior must be non-negative")
	
	p <- ncol(design)
	p0 <- p-1
	ngenes <- nrow(y)
	n <- ncol(y)
	d <- n-p

	if(length(contrast) == 1) {
		u <- rep.int(0,p)
		u[contrast] <- 1
		contrast <- u
	}
	if(length(contrast) != p) stop("length of contrast must match column dimension of design")
	if(all(contrast==0)) stop("contrast all zero")
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")

	if(!is.null(array.weights)) {
		if(any(array.weights <= 0)) stop("array.weights must be positive")
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		design <- design*sqrt(array.weights)
		y <- t(t(y)*sqrt(array.weights))
	}

	if(!is.null(block)) {
		block <- as.vector(block)
		if (length(block) != n) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,n,nblocks) == matrix(ub,n,nblocks,byrow = TRUE)
		cormatrix <- Z %*% (correlation * t(Z))
		diag(cormatrix) <- 1
		R <- chol(cormatrix)
		y <- t(backsolve(R, t(y), transpose = TRUE))
		design <- backsolve(R, design, transpose = TRUE)
 	}

	set.statistic <- match.arg(set.statistic,c("mean","floormean","mean50","msq"))

#	Reform design matrix so that contrast of interest is last column
	qr <- qr(contrast)
	Q <- qr.Q(qr,complete=TRUE)
	sign1 <- sign(qr$qr[1,1])
	Q <- cbind(Q[,-1],Q[,1])
	X <- design %*% Q
	qr <- qr(X)
	sign2 <- sign(qr$qr[p,p])
	signc <- sign1 * sign2

	if(is.null(var.prior) || is.null(df.prior)) {
#		Fit model to all genes
		effects <- qr.qty(qr,t(y))
#		Estimate global parameters s0 and d0
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(trend.var) covariate <- rowMeans(y,na.rm=TRUE) else covariate <- NULL
		sv <- squeezeVar(s2,df=d,covariate=covariate)
		d0 <- sv$df.prior
		s02 <- sv$var.prior
		if(trend.var) s02 <- s02[iset]
		effects <- effects[,iset,drop=FALSE]
		sd.post <- sqrt(sv$var.post[iset])
	} else {
		d0 <- df.prior
		s02 <- var.prior
		if(length(s02)>1) {
			names(s02) <- rownames(y)
			s02 <- s02[iset]
		}
		y <- y[iset,,drop=FALSE]
		effects <- qr.qty(qr,t(y))
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(is.finite(d0))
			sd.post <- sqrt( (d0*s02+d*s2)/(d0+d) )
		else
			sd.post <- sqrt(s02)
	}

#	From here, all results are for set only
	nset <- ncol(effects)
	if(p0>0)
		Y <- effects[-(1:p0),,drop=FALSE]
	else
		Y <- effects
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post

	statobs <- p <- rep(0,3)
	names(statobs) <- names(p) <- c("down","up","mixed")
	statrot <- array(0,c(nrot,3),dimnames=list(NULL,names(p)))

#	Convert to z-scores
	modt <- zscoreT(modt,df=d0+d)

#	Active proportions	
	if(!is.null(gene.weights)) {
		lgw <- length(gene.weights)
		if(lgw > nset && lgw==ngenes) {
			gene.weights <- gene.weights[iset]
		} else {
			if(lgw != nset) stop("length of gene.weights disagrees with size of set")
		}
		s <- sign(gene.weights)
		r1 <- mean(s*modt > sqrt(2))
		r2 <- mean(s*modt < -sqrt(2))
	} else {
		r1 <- mean(modt > sqrt(2))
		r2 <- mean(modt < -sqrt(2))
	}

#	Random rotations
	R <- matrix(rnorm(nrot*(d+1)),nrot,d+1)
	R <- R/sqrt(rowSums(R^2))
	Br <- R %*% Y
	s2r <- (matrix(YY,nrot,nset,byrow=TRUE)-Br^2)/d
	if(is.finite(d0))
		sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
	else
		sdr.post <- sqrt(s02)
	modtr <- signc*Br/sdr.post
	modtr <- zscoreT(modtr,df=d0+d)

	if(set.statistic=="msq") {
#		Observed statistics
		modt2 <- modt^2
		if(!is.null(gene.weights)) {
			modt2 <- abs(gene.weights)*modt2
			modt <- gene.weights*modt
		}
		statobs["mixed"] <- mean(modt2)
		statobs["up"] <- sum(modt2[modt > 0])/nset
		statobs["down"] <- sum(modt2[modt < 0])/nset
#		Simulated statistics   
		if(!is.null(gene.weights)) {
			gene.weights <- sqrt(abs(gene.weights))
			modtr <- t(gene.weights*t(modtr))
		}
		statrot[,"mixed"] <- rowMeans(modtr^2)
		statrot[,"up"] <- rowMeans(pmax(modtr,0)^2)
		statrot[,"down"] <- rowMeans(pmax(-modtr,0)^2)
#		p-values
		p <- (rowSums(t(statrot) >= statobs)+1)/(nrot+1)
	}

	if(set.statistic=="mean50") { 
#		Observed statistics
		if(!is.null(gene.weights)) modt <- gene.weights*modt
		half <- nset %/% 2L +1L
		meanTop <- function(x,n) mean(sort(x,partial=n)[n:length(x)])
		statobs["mixed"] <- meanTop(abs(modt),half)
		statobs["up"] <- meanTop(modt,half)
		statobs["down"] <- meanTop(-modt,half)
#		Simulated statistics
		if(!is.null(gene.weights)) modtr <- t(gene.weights*t(modtr))
		statrot[,"mixed"] <- apply(abs(modtr),1,meanTop,n=half)
		statrot[,"up"] <- apply(modtr,1,meanTop,n=half)
		statrot[,"down"] <- apply(-modtr,1,meanTop,n=half)
#		p-values
		p <- (rowSums(t(statrot) >= statobs) + 1)/(nrot + 1)
	}
   
	if(set.statistic=="floormean") { 
#		Observed statistics
		chimed <- qchisq(0.5,df=1)
		amodt <- pmax(abs(modt),chimed)
		if(!is.null(gene.weights)) {
			amodt <- gene.weights*amodt
			modt <- gene.weights*modt
		}
		statobs["mixed"] <- mean(amodt)
		statobs["up"] <- mean(pmax(modt,0))
		statobs["down"] <- mean(pmax(-modt,0))
#		Simulated statistics
		amodtr <- pmax(abs(modtr),chimed)
		if(!is.null(gene.weights)) {
			amodtr <- t(gene.weights*t(amodtr))
			modtr <- t(gene.weights*t(modtr))
		}
		statrot[,"mixed"] <- rowMeans(amodtr)
		statrot[,"up"] <- rowMeans(pmax(modtr,0))
		statrot[,"down"] <- rowMeans(pmax(-modtr,0))
#		p-values
		p <- (rowSums(t(statrot) >= statobs) + 1)/(nrot + 1)
	}
   
	if(set.statistic=="mean") { 
#		Observed statistics
		if(!is.null(gene.weights)) modt <- gene.weights*modt
		m <- mean(modt)
		statobs["mixed"] <- mean(abs(modt))
		statobs["up"] <- m
		statobs["down"] <- -m
#		Simulated statistics
		if(!is.null(gene.weights)) modtr <- t(gene.weights*t(modtr))
		m <- rowMeans(modtr)
		statrot[,"mixed"] <- rowMeans(abs(modtr))
		statrot[,"up"] <- m
		statrot[,"down"] <- -m
#		p-values
		p <- (rowSums(t(statrot) >= statobs) + 1)/(nrot + 1)
	}

#	Output
	out <- data.frame(c(r2,r1,r1+r2),p)
	dimnames(out) <- list(c("Down","Up","Mixed"),c("Active.Prop","P.Value"))
	new("Roast",list(p.value=out,var.prior=s02,df.prior=d0))
}


mroast <- function(iset=NULL,y,design,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,adjust.method="BH",midp=TRUE)
#  Rotation gene set testing with multiple sets
#  Gordon Smyth and Di Wu
#  Created 28 Jan 2010. Last revised 3 Feb 2012.
{ 
	if(is.null(iset)) iset <- rep(TRUE,nrow(y))
	if(!is.list(iset)) iset <- list(set = iset)
	nsets <- length(iset)
	if(is.null(names(iset))) names(iset) <- paste("set",1:nsets,sep="")
	if(!is.null(gene.weights)) if(length(gene.weights) != nrow(y)) stop("gene.weights must have length equal to nrow(y)")

#	Estimate var.prior and df.prior if not preset
	fit <- lmFit(y,design=design,weights=array.weights,block=block,correlation=correlation)
	covariate <- NULL
	if(trend.var) {
		covariate <- fit$Amean
		if(is.null(covariate)) covariate <- rowMeans(as.matrix(y))
	}
	sv <- squeezeVar(fit$sigma^2,df=fit$df.residual,covariate=covariate)
	var.prior <- sv$var.prior
	df.prior <- sv$df.prior

	pv <- adjpv <- active <- array(0,c(nsets,3),dimnames=list(names(iset),c("Down","Up","Mixed")))
	if(nsets<1) return(pv)
	for(i in 1:nsets) {
		out <- roast(iset=iset[[i]],y=y,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,nrot=nrot)[[1]]
		pv[i,] <- out$P.Value
		active[i,] <- out$Active.Prop
	}

#	Use mid-p-values or ordinary p-values?
	pv2 <- pv
	if(midp) pv2 <- pv2-1/2/(nrot+1)

	adjpv[,"Mixed"] <- p.adjust(pv2[,"Mixed"], method=adjust.method) 
	adjpv[,"Up"] <- p.adjust(pv2[,"Up"], method=adjust.method) 
	adjpv[,"Down"] <- p.adjust(pv2[,"Down"], method=adjust.method) 
	list(P.Value=pv, Adj.P.Value=adjpv, Active.Proportion=active) 
}

