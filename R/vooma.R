vooma <- function(y,design=NULL,correlation,block=NULL,plot=FALSE,span=NULL)
# Linear modelling of microarray data mean-variance modelling at the observational level.
# Creates an EList object for entry to lmFit() etc in the limma pipeline.
# Gordon Smyth and Charity Law
# Created 31 July 2012.  Last modified 5 Aug 2012.
{
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	narrays <- ncol(y)
	ngenes <- nrow(y)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		design <- matrix(1,narrays,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "GrandMean"
	}

#	Fit linear model
	if(is.null(block)) {
		fit <- lm.fit(design,t(y$E))
		mu <- fit$fitted.values
	} else {
		block <- as.vector(block)
		if(length(block)!=narrays) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,narrays,nblocks)==matrix(ub,narrays,nblocks,byrow=TRUE)
		cormatrix <- Z%*%(correlation*t(Z))
		diag(cormatrix) <- 1
		cholV <- chol(cormatrix)
		z <- backsolve(cholV,t(y$E),transpose=TRUE)
		X <- backsolve(cholV,design,transpose=TRUE)
		fit <- lm.fit(X,z)
		mu <- crossprod(cholV,fit$fitted.values)
	}
	s2 <- colMeans(fit$effects[-(1:fit$rank),,drop=FALSE]^2)

#	Fit lowess trend to sqrt-standard-deviations by ave log intensity
	sx <- rowMeans(y$E)
	sy <- sqrt(sqrt(s2))
	if(is.null(span)) if(ngenes<=10) span <- 1 else span <- 0.3+0.7*(10/ngenes)^0.5
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="Average log2 expression",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lines(l,col="red")
	}

#	Make interpolating rule
	f <- approxfun(l, rule=2)

#	Apply trend to individual observations
	w <- 1/f(t(mu))^4
	dim(w) <- dim(y)
	colnames(w) <- colnames(y)
	rownames(w) <- rownames(y)

#	Output
	y$meanvar.trend <- list(x=sx,y=sy)
	y$weights <- w
	y$design <- design
	y$span <- span
	y
}


voomaByGroup <- function(y,group,design=NULL,correlation,block=NULL,plot=FALSE,span=NULL,col=NULL)
#	Voom by group
#  Linear modelling of microarray data mean-variance modelling at the observational level by fitting group-specific trends.
#  Creates an EList object for entry to lmFit() etc in the limma pipeline.
#  Charity Law and Gordon Smyth
#  Created 13 February2013.
{
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))

#	Check group
	group <- as.factor(group)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		design <- model.matrix(~group)
		colnames(design) <- c("Int", levels(group)[-1])
	}

#	Check color
	if(is.null(col)) col <- rainbow(nlevels(group))
	color <- rep(col,length=nlevels(group))

	w <- matrix(nrow=nrow(y), ncol=ncol(y))
	colnames(w) <- colnames(y)
	rownames(w) <- rownames(y)
	sx <- matrix(nrow=nrow(y), ncol=nlevels(group))
	colnames(sx) <- levels(group)
	rownames(sx) <- rownames(y)
	sy <- matrix(nrow=nrow(y), ncol=nlevels(group))
	colnames(sy) <- levels(group)
	rownames(sy) <- rownames(y)
		for (lev in 1:nlevels(group)) {
		i <- group==levels(group)[lev]
		yi <- y[,i]
		designi <- design[i,,drop=FALSE]
		if(!is.null(block)) {
			blocki <- block[i]
			voomi <- vooma(y=yi,design=designi,correlation=correlation, block=blocki, plot=FALSE, span=span)
		} else {
			voomi <- vooma(y=yi,design=designi,plot=FALSE,span=span)
		}
		w[,i] <- voomi$weights
		sx[,lev] <- voomi$meanvar.trend$x
		sy[,lev] <- voomi$meanvar.trend$y	
	}
	span <- voomi$span

# 	Voom plot	
	if(plot) {		
		lev	<- 1
		colori <- color[lev]
		colori.transparent <- paste(rgb(col2rgb(colori)[1], col2rgb(colori)[2], col2rgb(colori)[3], maxColorValue=255), "10", sep="")
		plot(sx[,lev],sy[,lev],xlab="Average log2 expression",ylab="Sqrt( standard deviation )", xlim=c(min(as.vector(sx))*0.99, max(as.vector(sx))*1.01),ylim=c(min(as.vector(sy))*0.99, max(as.vector(sy))*1.01),pch=16, cex=0.25, col=colori.transparent)
		title("voom: Mean-variance trend")	
		l <- lowess(sx[,lev],sy[,lev],f=span)
		lines(l,col=colori)
		for (lev in 2:nlevels(group)) {
			colori <- color[lev]
			colori.transparent <- paste(rgb(col2rgb(colori)[1], col2rgb(colori)[2], col2rgb(colori)[3], maxColorValue=255), "10", sep="")
			points(sx[,lev], sy[,lev], 	pch=15, cex=0.25, col=colori.transparent)
			l <- lowess(sx[,lev],sy[,lev],f=span)
			lines(l,col=colori)	
		}
		legend("topright", levels(group), col=color, lty=1)
	}

# 	Output	
	y$meanvar.trend <- list(x=sx,y=sy)
	y$weights <- w
	y$design <- design
	y$span <- span
	y
}
