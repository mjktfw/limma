#  TOPTABLE.R

topTable <- function(fit,coef=NULL,number=10,genelist=fit$genes,adjust.method="BH",sort.by="B",resort.by=NULL,p.value=1,lfc=0,confint=FALSE)
#	Summary table of top genes, object-orientated version
#	Gordon Smyth
#	4 August 2003.  Last modified 19 Jan 2012.
{
	if(is.null(coef)) coef <- 1:ncol(fit)
	if(length(coef)>1) {
		coef <- unique(coef)
		if(length(fit$coef[1,coef]) < ncol(fit)) fit <- eBayes(fit[,coef])
		if(sort.by=="B") sort.by <- "F"
		return(topTableF(fit,number=number,genelist=genelist,adjust.method=adjust.method,sort.by=sort.by,p.value=p.value,lfc=lfc))
	}
	fit <- unclass(fit)
	ebcols <- c("t","p.value","lods")
	if(confint) ebcols <- c("s2.post","df.total",ebcols)
	toptable(fit=fit[c("coefficients","stdev.unscaled")],
		coef=coef,
		number=number,
		genelist=genelist,
		A=fit$Amean,
		eb=fit[ebcols],
		adjust.method=adjust.method,
		sort.by=sort.by,
		resort.by=resort.by,
		p.value=p.value,
		lfc=lfc,
		confint=confint)
}

topTableF <- function(fit,number=10,genelist=fit$genes,adjust.method="BH",sort.by="F",p.value=1,lfc=0)
#	Summary table of top genes by F-statistic
#	Gordon Smyth
#	27 August 2006. Last modified 19 Jan 2012.
{
#	Check input
	if(is.null(fit$F)) stop("F-statistics not available. Try topTable for individual coef instead.")
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ProbeID=genelist,stringsAsFactors=FALSE)
	sort.by <- match.arg(sort.by,c("F","none"))

	M <- as.matrix(fit$coefficients)
	if(is.null(colnames(M))) colnames(M) <- paste("Coef",1:ncol(M),sep="")
	rn <- rownames(M)
	if(is.null(rn) || length(rn) > length(unique(rn))) rn <- 1:nrow(M)
	adj.P.Value <- p.adjust(fit$F.p.value,method=adjust.method)

	if(lfc > 0) {
		big <- apply(abs(M)>lfc,1,any)
		big[is.na(big)] <- FALSE
		if(any(!big)) {
			fit <- fit[big,]
			if(!is.null(genelist)) genelist <- genelist[big,,drop=FALSE]
			M <- M[big,]
			rn <- rn[big]
			adj.P.Value <- adj.P.Value[big]
		}
	}
	if(nrow(fit) < number) number <- nrow(fit)

	if(sort.by=="F")
		o <- order(fit$F.p.value,decreasing=FALSE)[1:number]
	else
		o <- 1:number

	if(p.value < 1) {
		nsig <- sum(adj.P.Value[o] <= p.value,na.rm=TRUE)
		if(nsig < number) o <- o[1:nsig]
	}

	if(is.null(genelist))
		tab <- data.frame(M[o,,drop=FALSE])
	else
		tab <- data.frame(genelist[o,,drop=FALSE],M[o,,drop=FALSE])
	if(!is.null(fit$Amean)) tab <- data.frame(tab,AveExpr=fit$Amean[o])
	tab <- data.frame(tab,F=fit$F[o],P.Value=fit$F.p.value[o],adj.P.Val=adj.P.Value[o])
	rownames(tab) <- rn[o]
	tab
}

toptable <- function(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="BH",sort.by="B",resort.by=NULL,p.value=1,lfc=0,confint=FALSE,...)
#	Summary table of top genes
#	Gordon Smyth
#	21 Nov 2002. Last revised 30 Jan 2011.
{
#	Check input
	if(length(coef)>1) coef <- coef[1]
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ID=genelist,stringsAsFactors=FALSE)
	if(is.null(eb)) {
		fit$coefficients <- as.matrix(fit$coefficients)[,coef]
		fit$stdev.unscaled <- as.matrix(fit$stdev.unscaled)[,coef]
		eb <- ebayes(fit,...)
		coef <- 1
	}
	M <- as.matrix(fit$coefficients)[,coef]
	if(is.null(A)) {
		if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
	} else {
		if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
	}
	tstat <- as.matrix(eb$t)[,coef]
	P.Value <- as.matrix(eb$p.value)[,coef]
	B <- as.matrix(eb$lods)[,coef]
	rownum <- 1:length(M)
	sort.by <- match.arg(sort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B","none"))
	if(sort.by=="M") sort.by="logFC"
	if(sort.by=="A" || sort.by=="Amean") sort.by <- "AveExpr"
	if(sort.by=="T") sort.by <- "t"
	if(sort.by=="p") sort.by <- "P"

#	Apply multiple testing adjustment
	adj.P.Value <- p.adjust(P.Value,method=adjust.method)

#	Apply p.value or lfc thresholds	
	if(p.value < 1 | lfc > 0) {
		sig <- (adj.P.Value < p.value) & (abs(M) > lfc)
		if(any(is.na(sig))) sig[is.na(sig)] <- FALSE
		if(all(!sig)) return(data.frame())
		genelist <- genelist[sig,,drop=FALSE]
		M <- M[sig]
		A <- A[sig]
		tstat <- tstat[sig]
		P.Value <- P.Value[sig]
		adj.P.Value <- adj.P.Value[sig]
		B <- B[sig]
		rownum <- rownum[sig]
	}
	ord <- switch(sort.by,
		logFC=order(abs(M),decreasing=TRUE),
		AveExpr=order(A,decreasing=TRUE),
		P=order(P.Value,decreasing=FALSE),
		t=order(abs(tstat),decreasing=TRUE),
		B=order(B,decreasing=TRUE),
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
	if(confint) {
		margin.error <- sqrt(eb$s2.post[top])*fit$stdev.unscaled[top]*qnorm(0.975)
		tab$CI.025 <- M[top]-margin.error
		tab$CI.975 <- M[top]+margin.error
	}
	if(!is.null(A)) tab <- data.frame(tab,AveExpr=A[top])
	tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],adj.P.Val=adj.P.Value[top],B=B[top])
	rownames(tab) <- as.character(rownum)[top]
	if(!is.null(resort.by)) {
		resort.by <- match.arg(resort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B"))
		if(resort.by=="M") resort.by <- "logFC"
		if(resort.by=="A" || resort.by=="Amean") resort.by <- "AveExpr"
		if(resort.by=="p") resort.by <- "P"
		if(resort.by=="T") resort.by <- "t"
		ord <- switch(resort.by,
			logFC=order(tab$logFC,decreasing=TRUE),
			AveExpr=order(tab$AveExpr,decreasing=TRUE),
			P=order(tab$P.Value,decreasing=FALSE),
			t=order(tab$t,decreasing=TRUE),
			B=order(tab$B,decreasing=TRUE)
		)
		tab <- tab[ord,]
	}
#	attr(tab,"coef") <- coef
#	attr(tab,"adjust.method") <- adjust.method
	tab
}

