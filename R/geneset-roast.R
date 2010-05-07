##  ROAST.R

roast <- function(iset=NULL,y,design,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,nrot=999)
# rotation gene set testing for linear models
# Gordon Smyth and Di Wu
# 24 Apr 2008. Revised 7 May 2010.
{
	if(is.null(iset)) iset <- rep(TRUE,nrow(y))
	y <- as.matrix(y)
	design <- as.matrix(design)
	if(!is.null(df.prior) && df.prior<0) stop("df.prior must be non-negative")
	
	p <- ncol(design)
	p0 <- p-1
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
		sv <- squeezeVar(s2,df=d)
		d0 <- sv$df.prior
		s02 <- sv$var.prior
		effects <- effects[,iset,drop=FALSE]
		sd.post <- sqrt(sv$var.post[iset])
	} else {
		y <- y[iset,,drop=FALSE]
		effects <- qr.qty(qr,t(y))
		d0 <- df.prior
		s02 <- var.prior
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(is.finite(d0))
			sd.post <- sqrt( (d0*s02+d*s2)/(d0+d) )
		else
			sd.post <- sqrt(s02)
	}

#	From here, all results are for set only
	nset <- ncol(effects)
	Y <- effects[-(1:p0),,drop=FALSE]
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post

	statobs <- p <- rep(0,3)
	names(statobs) <- names(p) <- c("mixed","up","down")
	statrot <- array(0,c(nrot,3),dimnames=list(NULL,names(p)))

#	Convert to z-scores
	modt <- zscoreT(modt,df=d0+d)

#	Active proportions	
	if(!is.null(gene.weights)) {
		if(length(gene.weights) != nset) stop("length of gene.weights disagrees with size of set")
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
		if(!is.null(gene.weights)) modt2 <- gene.weights*modt2
		statobs["mixed"] <- mean(modt2)
		statobs["up"] <- sum(modt2[modt > 0])/nset
		statobs["down"] <- sum(modt2[modt < 0])/nset
#		Simulated statistics   
		if(!is.null(gene.weights)) {
			gene.weights <- s*sqrt(abs(gene.weights))
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
	out <- data.frame(ActiveProp=c(r1+r2,r1,r2),P.Value=p)
	row.names(out) <- c("mixed","up","down")
	out
}

