diffSplice <- function(fit,geneid,exonid=NULL,robust=FALSE,verbose=TRUE)
#	Test for splicing variants between conditions
#	using linear model fit of exon data.
#	Gordon Smyth and Charity Law
#	Created 13 Dec 2013.  Last modified 22 Sep 2014.
{
	exon.genes <- fit$genes
	if(is.null(exon.genes)) exon.genes <- data.frame(ExonID=1:nrow(fit))
	if(length(geneid)==1) {
		genecolname <- as.character(geneid)
		geneid <- exon.genes[[genecolname]]
	} else {
		exon.genes$GeneID <- geneid
		genecolname <- "GeneID"
	}
	if(!is.null(exonid))
		if(length(exonid)==1) {
			exoncolname <- as.character(exonid)
			exonid <- exon.genes[[exoncolname]]
		} else {
			exon.genes$ExonID <- exonid
			exoncolname <- "ExonID"
		}
	else
		exoncolname <- NULL

#	Sort by geneid
	if(is.null(exonid))
		o <- order(geneid)
	else
		o <- order(geneid,exonid)
	geneid <- geneid[o]
	exon.genes <- exon.genes[o,,drop=FALSE]
	exon.coefficients <- fit$coefficients[o,,drop=FALSE]
	exon.stdev.unscaled <- fit$stdev.unscaled[o,,drop=FALSE]
	exon.df.residual <- fit$df.residual[o]
	exon.s2 <- fit$sigma[o]^2

# 	Count exons and get genewise variances
	exon.stat <- cbind(1,exon.df.residual,exon.s2)
	gene.sum <- rowsum(exon.stat,geneid,reorder=FALSE)
	gene.nexons <- gene.sum[,1]
	gene.df.residual <- gene.sum[,2]
	gene.s2 <- gene.sum[,3] / gene.sum[,1]
	if(verbose) {
		cat("Total number of exons: ", length(geneid), "\n")
		cat("Total number of genes: ", length(gene.nexons), "\n")
		cat("Number of genes with 1 exon: ", sum(gene.nexons==1), "\n")
		cat("Mean number of exons in a gene: ", round(mean(gene.nexons),0), "\n")
		cat("Max number of exons in a gene: ", max(gene.nexons), "\n")
	}

#	Posterior genewise variances
	squeeze <- squeezeVar(var=gene.s2, df=gene.df.residual, robust=robust)

#	Remove genes with only 1 exon
	gene.keep <- gene.nexons>1
	ngenes <- sum(gene.keep)
	if(ngenes==0) stop("No genes with more than one exon")

	exon.keep <- rep(gene.keep,gene.nexons)
	geneid <- geneid[exon.keep]
	exon.genes <- exon.genes[exon.keep,,drop=FALSE]
	exon.coefficients <- exon.coefficients[exon.keep,,drop=FALSE]
	exon.stdev.unscaled <- exon.stdev.unscaled[exon.keep,,drop=FALSE]
	exon.df.residual <- exon.df.residual[exon.keep]

	gene.nexons <- gene.nexons[gene.keep]
	gene.df.test <- gene.nexons-1
	gene.df.residual <- gene.df.residual[gene.keep]
	if(robust) squeeze$df.prior <- squeeze$df.prior[gene.keep]
	gene.df.total <- gene.df.residual+squeeze$df.prior
	gene.df.total <- pmin(gene.df.total,sum(gene.df.residual))
	gene.s2.post <- squeeze$var.post[gene.keep]

# 	Genewise betas
	u2 <- 1/exon.stdev.unscaled^2
	u2.rowsum <- rowsum(u2,geneid,reorder=FALSE)
	gene.betabar <- rowsum(exon.coefficients*u2,geneid,reorder=FALSE) / u2.rowsum

#	T-statistics for exon-level tests
	g <- rep(1:ngenes,times=gene.nexons)
	exon.coefficients <- exon.coefficients-gene.betabar[g,,drop=FALSE]
	exon.t <- exon.coefficients / exon.stdev.unscaled / sqrt(gene.s2.post[g])
	gene.F <- rowsum(exon.t^2,geneid,reorder=FALSE) / gene.df.test
	exon.1mleverage <- 1 - (u2 / u2.rowsum[g,,drop=FALSE])
	exon.coefficients <- exon.coefficients / exon.1mleverage
	exon.t <- exon.t / sqrt(exon.1mleverage)
	exon.p.value <- 2 * pt(abs(exon.t), df=gene.df.total[g], lower.tail=FALSE)
	gene.F.p.value <- pf(gene.F, df1=gene.df.test, df2=gene.df.total, lower.tail=FALSE)

#	Exon level output
	out <- new("MArrayLM",list())
	out$genes <- exon.genes
	out$genecolname <- genecolname
	out$exoncolname <- exoncolname
	out$coefficients <- exon.coefficients
	out$t <- exon.t
	out$p.value <- exon.p.value

#	Gene level output
	out$gene.df.prior <- squeeze$df.prior
	out$gene.df.residual <- gene.df.residual
	out$gene.df.total <- gene.df.total
	out$gene.s2.post <- gene.s2.post
	out$gene.F <- gene.F
	out$gene.F.p.value <- gene.F.p.value

#	Which columns of exon.genes contain gene level annotation?
	gene.lastexon <- cumsum(gene.nexons)
	gene.firstexon <- gene.lastexon-gene.nexons+1
	no <- logical(nrow(exon.genes))
	isdup <- vapply(exon.genes,duplicated,no)[-gene.firstexon,,drop=FALSE]
	isgenelevel <- apply(isdup,2,all)
	out$gene.genes <- exon.genes[gene.lastexon,isgenelevel, drop=FALSE]
	out$gene.genes$NExons <- gene.nexons
	out$gene.firstexon <- gene.firstexon
	out$gene.lastexon <- gene.lastexon

#	Simes adjustment of exon level p-values
	penalty <- rep_len(1L,length(g))
	penalty[gene.lastexon] <- 1L-gene.nexons
	penalty <- cumsum(penalty)[-gene.lastexon]
	penalty <- penalty / rep(gene.nexons-1L,gene.nexons-1L)
	g2 <- g[-gene.lastexon]

	out$gene.simes.p.value <- gene.F.p.value
	for (j in 1:ncol(fit)) {
		o <- order(g,exon.p.value[,j])
		p.adj <- pmin(exon.p.value[o,j][-gene.lastexon] / penalty, 1)
		o <- order(g2,p.adj)
		out$gene.simes.p.value[,j] <- p.adj[o][gene.firstexon-0L:(ngenes-1L)]
	}

	out
}

topSplice <- function(fit, coef=ncol(fit), test="simes", number=10, FDR=1)
#	Collate diffSplice results into data.frame, ordered from most significant at top
#	Gordon Smyth
#	Created 18 Dec 2013.  Last modified 18 Aug 2014.
{
	coef <- coef[1]
	test <- match.arg(test,c("simes","F","f","t"))
	if(test=="f") test <- "F"
	switch(test,
	"t" = {
		number <- min(number,nrow(fit$coefficients))
		P <- fit$p.value[,coef]
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number,sum(BH<FDR))
		o <- order(P)[1:number]
		data.frame(fit$genes[o,,drop=FALSE],logFC=fit$coefficients[o,coef],t=fit$t[o,coef],P.Value=P[o],FDR=BH[o])
	},
	F = {
		number <- min(number,nrow(fit$gene.F))
		P <- fit$gene.F.p.value[,coef]
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number,sum(BH<FDR))
		o <- order(P)[1:number]
		data.frame(fit$gene.genes[o,,drop=FALSE],F=fit$gene.F[o,coef],P.Value=P[o],FDR=BH[o])
	},
	simes = {
		number <- min(number,nrow(fit$gene.F))
		P <- fit$gene.simes.p.value[,coef]
		BH <- p.adjust(P, method="BH")
		if(FDR<1) number <- min(number,sum(BH<FDR))
		o <- order(P)[1:number]
		data.frame(fit$gene.genes[o,,drop=FALSE],P.Value=P[o],FDR=BH[o])
	}
	)
}

plotSplice <- function(fit, coef=ncol(fit), geneid=NULL, genecolname=NULL, rank=1L, FDR = 0.05)
#	Plot exons of chosen gene
#	fit is output from diffSplice
#	Gordon Smyth and Yifang Hu
#	Created 3 Jan 2014.  Last modified 7 October 2014.
{
	if(is.null(genecolname)) 
		genecolname <- fit$genecolname
	else
		genecolname <- as.character(genecolname)

	if(is.null(geneid)) {
#		Find gene from specified rank 
		if(rank==1L)
			i <- which.min(fit$gene.F.p.value[,coef])
		else
			i <- order(fit$gene.F.p.value[,coef])[rank]
		geneid <- paste(fit$gene.genes[i,genecolname], collapse = ".")
	} else {
#		Find gene from specified name
		geneid <- as.character(geneid)
		i <- which(fit$gene.genes[,genecolname]==geneid)[1]
		if(!length(i)) stop(paste("geneid",geneid,"not found"))
	}

#	Row numbers containing exons
	j <- fit$gene.firstexon[i]:fit$gene.lastexon[i]

	exoncolname <- fit$exoncolname

#	Get strand if possible		
	strcol <- grepl("strand", colnames(fit$gene.genes), ignore.case = TRUE)
	if(any(strcol)) geneid <- paste0(geneid, " (", as.character(fit$gene.genes[i, strcol])[1], ")")

	if(is.null(exoncolname)) {

		plot(fit$coefficients[j,coef], xlab = "Exon", ylab = "logFC (this exon vs rest)", main = geneid, type = "b")

	} else {

		exon.id <- fit$genes[j, exoncolname]
		xlab <- paste("Exon", exoncolname, sep = " ")

		plot(fit$coefficients[j, coef], xlab = "", ylab = "logFC (this exon vs rest)", main = geneid,type = "b", xaxt = "n")
		axis(1, at = 1:length(j), labels = exon.id, las = 2, cex.axis = 0.5)
		mtext(xlab, side = 1, padj = 5.2)

#		Mark the topSpliced exons
		top <- topSplice(fit, coef = coef, number = Inf, test = "t", FDR = FDR)
		m <- which(top[,genecolname] %in% fit$gene.genes[i,genecolname])
		if(length(m) > 0){
			if(length(m) == 1)
				cex <- 1.5
			else {
				abs.fdr <- abs(log10(top$FDR[m]))
				from <- range(abs.fdr)
				to <- c(1,2)
				cex <- (abs.fdr - from[1])/diff(from) * diff(to) + to[1]
			}	
			mark <- match(top[m, exoncolname], exon.id)
			points((1:length(j))[mark], fit$coefficients[j[mark], coef], col = "red", pch = 16, cex = cex)
		}

	}

	abline(h=0,lty=2)
	invisible()
}
