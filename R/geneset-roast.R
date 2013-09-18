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

roast <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999)
UseMethod("roast")

roast.EList <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999)
# Gordon Smyth
# Created 4 Jan 2013.  Last modified 22 Jan 2013.
{
	if(is.null(design)) design <- y$design
	if(is.null(weights) && is.null(array.weights)) weights <- y$weights
	y <- as.matrix(y)
	roast(y=y,index=index,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,trend.var=trend.var,nrot=nrot)
}

roast.MAList <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999)
# Gordon Smyth
# Created 4 Jan 2013.  Last modified 22 Jan 2013.
{
	if(is.null(design)) design <- y$design
	if(is.null(weights) && is.null(array.weights)) weights <- y$weights
	y <- as.matrix(y)
	roast(y=y,index=index,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,trend.var=trend.var,nrot=nrot)
}

roast.default <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999)
# Rotation gene set testing for linear models
# Gordon Smyth and Di Wu
# Created 24 Apr 2008.  Last modified 19 Sep 2013.
{
#	Check y
	y <- as.matrix(y)
	ngenes <- nrow(y)
	n <- ncol(y)

#	Check index
	if(is.list(index)) return(mroast(y=y,index=index,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,trend.var=trend.var,nrot=nrot))
	if(is.null(index)) index <- rep.int(TRUE,ngenes)

#	Check design
	if(is.null(design)) stop("no design matrix")
	design <- as.matrix(design)
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	p <- ncol(design)
	p0 <- p-1
	d <- n-p

#	Check contrast
	if(length(contrast) == 1) {
		u <- rep.int(0,p)
		u[contrast] <- 1
		contrast <- u
	}
	if(length(contrast) != p) stop("length of contrast must match column dimension of design")
	if(all(contrast==0)) stop("contrast all zero")

#	Check set.statistic
	set.statistic <- match.arg(set.statistic,c("mean","floormean","mean50","msq"))

#	Check var.prior and df.prior
	if(!is.null(var.prior) && var.prior<0) stop("var.prior must be non-negative")
	if(!is.null(df.prior) && df.prior<0) stop("df.prior must be non-negative")
	
#	Check and divide out array weights
	if(!is.null(array.weights)) {
		if(!is.null(weights)) stop("Can't specify both array.weights and observational weights")
		if(any(array.weights <= 0)) stop("array.weights must be positive")
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		design <- design*sqrt(array.weights)
		y <- t(t(y)*sqrt(array.weights))
	}

#	Check and divide out block correlation
	if(!is.null(block)) {
		if(!is.null(weights)) stop("Can't use block with weights")
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

#	Check weights
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		dimw <- dim(weights)
		if(dimw[1]!=ngenes || dimw[2]!=n) stop("weights must have same dimensions as y")
		if(any(weights<=0)) stop("weights must be positive")
	}

#	Reform design matrix so that contrast of interest is last column
#	qr <- qr(contrast)
#	Q <- qr.Q(qr,complete=TRUE)
#	sign1 <- sign(qr$qr[1,1])
#	Q <- cbind(Q[,-1],Q[,1])
#	X <- design %*% Q
	X <- contrastAsCoef(design, contrast, first=FALSE)$design
	qr <- qr(X)
	signc <- sign(qr$qr[p,p])

	if(is.null(var.prior) || is.null(df.prior)) {
#		Fit model to all genes
		if(is.null(weights)) {
			effects <- qr.qty(qr,t(y))
		} else {
			ws <- sqrt(weights)
			effects <- matrix(0,n,ngenes)
			signc <- rep.int(0,ngenes)
			for (g in 1:ngenes) {
				wX <- X*ws[g,]
				wy <- y[g,]*ws[g,]
				qrX <- qr(wX)
				signc[g] <- sign(qrX$qr[p,p])
				effects[,g] <- qr.qty(qrX,wy)
			}
			signc <- signc[index]
		}
#		Estimate global parameters s0 and d0
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(trend.var) covariate <- rowMeans(y,na.rm=TRUE) else covariate <- NULL
		sv <- squeezeVar(s2,df=d,covariate=covariate)
		d0 <- sv$df.prior
		s02 <- sv$var.prior
		if(trend.var) s02 <- s02[index]
		effects <- effects[,index,drop=FALSE]
		sd.post <- sqrt(sv$var.post[index])
	} else {
		d0 <- df.prior
		s02 <- var.prior
		if(length(s02)>1) {
			names(s02) <- rownames(y)
			s02 <- s02[index]
		}
		y <- y[index,,drop=FALSE]
		if(is.null(weights)) {
			effects <- qr.qty(qr,t(y))
		} else {
			ws <- sqrt(weights[index,,drop=FALSE])
			nset <- nrow(y)
			effects <- matrix(0,n,nset)
			signc <- rep.int(0,nset)
			for (g in 1:nset) {
				wX <- X*ws[g,]
				wy <- y[g,]*ws[g,]
				qrX <- qr(wX)
				signc[g] <- sign(qrX$qr[p,p])
				effects[,g] <- qr.qty(qrX,wy)
			}
		}
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
			gene.weights <- gene.weights[index]
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
	new("Roast",list(p.value=out,var.prior=s02,df.prior=d0,ngenes.in.set=nset))
}

mroast <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,adjust.method="BH",midp=TRUE,sort="directional")
UseMethod("mroast")

mroast.EList <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,adjust.method="BH",midp=TRUE,sort="directional")
# Gordon Smyth
# Created 8 Jan 2013.
{
	if(is.null(design)) design <- y$design
	if(is.null(weights) && is.null(array.weights)) weights <- y$weights
	y <- as.matrix(y)
	mroast(y=y,index=index,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,trend.var=trend.var,nrot=nrot)
}

mroast.MAList <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,adjust.method="BH",midp=TRUE,sort="directional")
# Gordon Smyth
# Created 8 Jan 2013.
{
	if(is.null(design)) design <- y$design
	if(is.null(weights) && is.null(array.weights)) weights <- y$weights
	y <- as.matrix(y)
	mroast(y=y,index=index,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,trend.var=trend.var,nrot=nrot)
}

mroast.default <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,adjust.method="BH",midp=TRUE,sort="directional")
#  Rotation gene set testing with multiple sets
#  Gordon Smyth and Di Wu
#  Created 28 Jan 2010. Last revised 19 Sep 2013.
{
#	Check y
	y <- as.matrix(y)

#	Check index
	if(is.null(index)) index <- rep(TRUE,nrow(y))
	if(!is.list(index)) index <- list(set = index)
	nsets <- length(index)
	if(is.null(names(index))) names(index) <- paste("set",1:nsets,sep="")

#	Check design
	if(is.null(design)) stop("No design matrix")
	design <- as.matrix(design)

#	Check gene.weights
	if(!is.null(gene.weights)) if(length(gene.weights) != nrow(y)) stop("gene.weights must have length equal to nrow(y)")

#	Check weights
	if(!is.null(weights)) {
		if(!is.null(array.weights)) stop("Can't specify both array weights and observational weights")
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(y))) stop("weights must have same dimensions as y")
	}

#	Check array.weights
	if(!is.null(array.weights)) {
		if(length(array.weights) != ncol(y)) stop("array.weights wrong length")
		weights <- array.weights
	}

#	Estimate var.prior and df.prior if not preset
	if(is.null(var.prior) || is.null(df.prior)) {
		fit <- lmFit(y,design=design,weights=weights,block=block,correlation=correlation)
		if(trend.var) {
			covariate <- fit$Amean
			if(is.null(covariate)) covariate <- rowMeans(y)
		} else {
			covariate=NULL
		}
		sv <- squeezeVar(fit$sigma^2,df=fit$df.residual,covariate=covariate)
		var.prior <- sv$var.prior
		df.prior <- sv$df.prior
	}

	pv <- adjpv <- active <- array(0,c(nsets,3),dimnames=list(names(index),c("Down","Up","Mixed")))
	NGenes <- rep(0,nsets)
	if(nsets<1) return(pv)
	for(i in 1:nsets) {
		out <- roast(y=y,index=index[[i]],design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,nrot=nrot)
		pv[i,] <- out$p.value$P.Value
		active[i,] <- out$p.value$Active.Prop
		NGenes[i] <- out$ngenes.in.set
	}

#	Use mid-p-values or ordinary p-values?
	pv2 <- pv
	if(midp) pv2 <- pv2-1/2/(nrot+1)

#	adjpv[,"Down"] <- p.adjust(pv2[,"Down"], method=adjust.method)
#	adjpv[,"Up"] <- p.adjust(pv2[,"Up"], method=adjust.method)
#	adjpv[,"Mixed"] <- p.adjust(pv2[,"Mixed"], method=adjust.method)
#	list(P.Value=pv, Adj.P.Value=adjpv, Active.Proportion=active)

#	New-style output
	Up <- pv[,"Up"] < pv[,"Down"]
	Direction <- rep.int("Down",nsets); Direction[Up] <- "Up"
	TwoSidedP <- pv[,"Down"]; TwoSidedP[Up] <- pv[Up,"Up"]; TwoSidedP <- 2*TwoSidedP
	TwoSidedP2 <- pv2[,"Down"]; TwoSidedP2[Up] <- pv2[Up,"Up"]; TwoSidedP2 <- 2*TwoSidedP2
	tab <- data.frame(
		NGenes=NGenes,
		PropDown=active[,"Down"],
		PropUp=active[,"Up"],
		Direction=Direction,
		PValue=TwoSidedP,
		FDR=p.adjust(TwoSidedP2,method="BH"),
		PValue.Mixed=pv[,"Mixed"],
		FDR.Mixed=p.adjust(pv2[,"Mixed"],method="BH"),
		row.names=names(index),
		stringsAsFactors=FALSE
	)

#	Sort by p-value
	sort <- match.arg(sort,c("directional","mixed","none"))
	if(sort=="none") return(tab)
	if(sort=="directional")
		o <- order(tab$PValue)
	else
		o <- order(tab$PValue.Mixed)
	tab[o,,drop=FALSE]
}

