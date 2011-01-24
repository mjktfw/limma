#  Plot regularized linear discriminant functions

plotRLDF <- function(y,design=NULL,z=NULL,labels.y=NULL,labels.z=NULL,col.y=1,col.z=1,df.prior=5,show.dimensions=c(1,2),main=NULL,nprobes=500,...)
#	Regularized linear discriminant function
#	Di Wu and Gordon Smyth
#	29 June 2009.  Last revised 30 June 2009.
{
	y <- as.matrix(y)
	g <- nrow(y)
	n <- ncol(y)
	if(is.null(design)) {
		if(is.null(labels.y)) stop("groups not specified")
		f <- as.factor(labels.y)
		design <- model.matrix(~f)
	}
	design <- as.matrix(design)
	if(nrow(design) != n) stop("nrow(design) doesn't match ncol(y)")
	if(is.null(labels.y)) labels.y <- 1:n

#	Project onto between and within spaces
#	Discard first column as intercept
	qrd <- qr(design)
	p <- qrd$rank
	if(p==n) stop("No residual degrees of freedom")
	U <- qr.qty(qrd, t(y))
	UB <- U[2:p,,drop=FALSE]
	UW <- U[(p+1):n,,drop=FALSE]

#	Prior variance
	s <- colMeans(UW*UW)
	s0 <- median(s)

#	Select probes by moderated F
	if(g>nprobes) {
		modF <- colMeans(UB*UB)/(s+df.prior*s0)
		o <- order(modF,decreasing=TRUE)
		top <- o[1:nprobes]
		y <- y[top,,drop=FALSE]
		if(!is.null(z)) z <- z[top,,drop=FALSE]
		UB <- UB[,top,drop=FALSE]
		UW <- UW[,top,drop=FALSE]
		g <- nprobes
	}

#	Within group SS
	W <- crossprod(UW)

#	Regularized within group SS
	Wreg <- W+diag(df.prior*s0,g,g)

#	Ratio of between to within SS
	WintoB <- backsolve(chol(Wreg),t(UB),transpose=TRUE)

#	Linear discriminant gene weights
	d1 <- show.dimensions[1]
	d2 <- show.dimensions[2]
	metagenes <- svd(WintoB,nu=max(d1,d2),nv=0)$u

	d1.y <- t(y)%*%metagenes[,d1]
	d2.y <- t(y)%*%metagenes[,d2]
	if(!is.null(z)) {
		z <- as.matrix(z)
		d1.z <- t(z)%*%metagenes[,d1]
		d2.z <- t(z)%*%metagenes[,d2]
		if(is.null(labels.z)) labels.z <- letters[1:ncol(z)]
	} else {
		d1.z <- d2.z <- numeric(0)
		labels.z <- character(0)
	}
	lab.y <- as.character(labels.y)
	lab.z <- as.character(labels.z)
	plot(c(d1.y,d1.z),c(d2.y,d2.z),type="n",main=main,
		xlab=paste("Discriminant Function",d1),
		ylab=paste("Discriminant Function",d2))
	text(d1.y,d2.y,labels=lab.y,col=col.y,...)
	if(!is.null(z)) text(d1.z,d2.z,labels=lab.z,col=col.z,...)
	invisible(list(training=cbind(d1=d1.y,d2=d2.y),predicting=cbind(d1=d1.z,d2=d2.z),metagenes=metagenes,top=top))
}

