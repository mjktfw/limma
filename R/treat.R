###  treat.R

treat <- function(fit, lfc=0)
#  Moderated t-statistics with threshold
#  Davis McCarthy, Gordon Smyth
#  25 July 2008. Last revised 3 November 2010.
{
	coefficients <- as.matrix(fit$coefficients)
	stdev.unscaled <- as.matrix(fit$stdev.unscaled)
	sigma <- fit$sigma
	df.residual <- fit$df.residual
	if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
		is.null(df.residual)) 
		stop("No data, or argument is not a valid lmFit object")
	if (all(df.residual == 0)) 
		stop("No residual degrees of freedom in linear model fits")
	if (all(!is.finite(sigma))) 
		stop("No finite residual standard deviations")
	sv <- squeezeVar(sigma^2, df.residual)
	fit$df.prior <- sv$df.prior
	fit$s2.prior <- sv$var.prior
	fit$s2.post <- sv$var.post
	df.total <- df.residual + sv$df.prior
	df.pooled <- sum(df.residual,na.rm=TRUE)
	df.total <- pmin(df.total,df.pooled)
	fit$df.total <- df.total
	lfc <- abs(lfc)
	acoef <- abs(coefficients)
	se <- stdev.unscaled*sqrt(fit$s2.post)
	tstat.right <- (acoef-lfc)/se
	tstat.left <- (acoef+lfc)/se
	fit$t <- array(0,dim(coefficients),dimnames=dimnames(coefficients))
	fit$p.value <- pt(tstat.right, df = df.total,lower=FALSE) + pt(tstat.left,df = df.total,lower=FALSE)
	tstat.right <- pmax(tstat.right,0)
	fc.up <- (coefficients > lfc)
	fc.down <- (coefficients < -lfc)
	fit$t[fc.up] <- tstat.right[fc.up]
	fit$t[fc.down] <- -tstat.right[fc.down]
	fit
}

topTreat <- function(fit,coef=1,number=10,genelist=fit$genes,adjust.method="BH",sort.by="p",resort.by=NULL,p.value=1)
#	Summary table of top genes by treat
#	Gordon Smyth
#	15 June 2009.  Last modified 17 March 2010.
{
#	Check input
	if(length(coef)>1) coef <- coef[1]
	M <- as.matrix(fit$coefficients)[,coef]
	A <- fit$Amean
	if(is.null(A)) {
		if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
	} else {
		if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
	}
	tstat <- as.matrix(fit$t)[,coef]
	P.Value <- as.matrix(fit$p.value)[,coef]
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ID=genelist,stringsAsFactors=FALSE)
	sort.by <- match.arg(sort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","none"))
	if(sort.by=="M") sort.by="logFC"
	if(sort.by=="A" || sort.by=="Amean") sort.by <- "AveExpr"
	if(sort.by=="T") sort.by <- "t"
	if(sort.by=="p") sort.by <- "P"

#	Apply multiple testing adjustment
	adj.P.Value <- p.adjust(P.Value,method=adjust.method)

#	Apply p.value threshold
	if(p.value < 1) {
		sig <- (adj.P.Value < p.value)
		if(any(is.na(sig))) sig[is.na(sig)] <- FALSE
		if(all(!sig)) return(data.frame())
		genelist <- genelist[sig,,drop=FALSE]
		M <- M[sig]
		A <- A[sig]
		tstat <- tstat[sig]
		P.Value <- P.Value[sig]
		adj.P.Value <- adj.P.Value[sig]
	}
	ord <- switch(sort.by,
		logFC=order(abs(M),decreasing=TRUE),
		AveExpr=order(A,decreasing=TRUE),
		P=order(P.Value,decreasing=FALSE),
		t=order(abs(tstat),decreasing=TRUE),
		none=1:length(M)
	)
	if(length(M) < number) number <- length(M)
	if(number < 1) return(data.frame())
	top <- ord[1:number]
	if(is.null(genelist))
		tab <- data.frame(logFC=M[top])
	else {
		tab <- data.frame(genelist[top,,drop=FALSE],logFC=M[top],stringsAsFactors=FALSE)
	}
	if(!is.null(A)) tab <- data.frame(tab,AveExpr=A[top])
	tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],adj.P.Val=adj.P.Value[top])
	rownames(tab) <- as.character(1:length(M))[top]
	if(!is.null(resort.by)) {
		resort.by <- match.arg(resort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t"))
		if(resort.by=="M") resort.by <- "logFC"
		if(resort.by=="A" || resort.by=="Amean") resort.by <- "AveExpr"
		if(resort.by=="p") resort.by <- "P"
		if(resort.by=="T") resort.by <- "t"
		ord <- switch(resort.by,
			logFC=order(tab$logFC,decreasing=TRUE),
			AveExpr=order(tab$AveExpr,decreasing=TRUE),
			P=order(tab$P.Value,decreasing=FALSE),
			t=order(tab$t,decreasing=TRUE)
		)
		tab <- tab[ord,]
	}
	tab
}

