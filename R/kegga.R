kegga <- function(de,...) UseMethod("kegga")

kegga.MArrayLM <- function(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, ...)
#	KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway analysis of DE genes from linear model fit
#	Gordon Smyth and Yifang Hu
#	Created 15 May 2015. Last modified 3 June 2015.
{
#	Avoid argument collision with default method
	dots <- names(list(...))
	if("universe" %in% dots) stop("kegga.MArrayLM defines its own universe",call.=FALSE)
	if((!is.logical(trend) || trend) && "covariate" %in% dots) stop("kegga.MArrayLM defines it own covariate",call.=FALSE)

#	Check fit
	if(is.null(de$coefficients)) stop("de does not appear to be a valid MArrayLM fit object.")
	if(is.null(de$p.value)) stop("p.value not found in fit object, perhaps need to run eBayes first.")	
	if(length(coef) != 1) stop("Only one coef can be specified.")
	ngenes <- nrow(de)

#	Check geneid
#	Can be either a vector of gene IDs or an annotation column name
	geneid <- as.character(geneid)
	if(length(geneid) == ngenes) {
		universe <- geneid
	} else {
		if(length(geneid) == 1L) {
			universe <- de$genes[[geneid]]
			if(is.null(universe)) stop("Column ",geneid," not found in de$genes")
		} else
			stop("geneid of incorrect length")
	}

#	Check trend
#	Can be logical, or a numeric vector of covariate values, or the name of the column containing the covariate values
	if(is.logical(trend)) {
		if(trend) {
			covariate <- de$Amean
			if(is.null(covariate)) stop("Amean not found in fit")
		}
	} else {
		if(is.numeric(trend)) {
			if(length(trend) != ngenes) stop("If trend is numeric, then length must equal nrow(de)")
			covariate <- trend
			trend <- TRUE
		} else {
			if(is.character(trend)) {
				if(length(trend) != 1L) stop("If trend is character, then length must be 1")
				covariate <- de$genes[[trend]]
				if(is.null(covariate)) stop("Column ",trend," not found in de$genes")
				trend <- TRUE
			} else
				stop("trend is neither logical, numeric nor character")
		}
	}

#	Check FDR
	if(!is.numeric(FDR) | length(FDR) != 1) stop("FDR must be numeric and of length 1.")
	if(FDR < 0 | FDR > 1) stop("FDR should be between 0 and 1.")

#	Get up and down DE genes
	fdr.coef <- p.adjust(de$p.value[,coef], method = "BH")
	EG.DE.UP <- universe[fdr.coef < FDR & de$coef[,coef] > 0]
	EG.DE.DN <- universe[fdr.coef < FDR & de$coef[,coef] < 0]
	DEGenes <- list(Up=EG.DE.UP, Down=EG.DE.DN)

#	If no DE genes, return data.frame with 0 rows
	if(length(EG.DE.UP)==0 && length(EG.DE.DN)==0) {
		message("No DE genes")
		return(data.frame())
	}

	if(trend)
		kegga(de=DEGenes, universe = universe, covariate=covariate, ...)
	else
		kegga(de=DEGenes, universe = universe, ...)
}

kegga.default <- function(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL, plot=FALSE, ...)
#	KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway analysis of DE genes
#	Gordon Smyth and Yifang Hu
#	Created 18 May 2015.  Modified 4 August 2015.
{
#	Ensure de is a list
	if(!is.list(de)) de <- list(DE = de)

#	Stop if components of de are not vectors
	if(!all(vapply(de,is.vector,TRUE))) stop("components of de should be vectors")

#	Ensure gene IDs are of character mode
	de <- lapply(de, as.character)
	if(!is.null(universe)) universe <- as.character(universe)

#	Ensure all gene sets have unique names
	nsets <- length(de)
	names(de) <- trimWhiteSpace(names(de))
	NAME <- names(de)
	i <- which(NAME == "" | is.na(NAME))
	NAME[i] <- paste0("DE",i)
	names(de) <- makeUnique(NAME)

#	Select species
	species <- match.arg(species, c("Hs", "Mm", "Rn", "Dm"))

#	Fit trend in DE with respect to the covariate, combining all de lists
	if(!is.null(covariate)) {
		covariate <- as.numeric(covariate)
		if(length(covariate) != length(covariate)) stop("universe and covariate must have same length")
		isDE <- as.numeric(universe %in% unlist(de))
		o <- order(covariate)
		prior.prob <- covariate
		span <- approx(x=c(20,200),y=c(1,0.5),xout=sum(isDE),rule=2)$y
		prior.prob[o] <- tricubeMovingAverage(isDE[o],span=span,full.length=TRUE)
		if(plot) barcodeplot(covariate, index=(isDE==1), worm=TRUE, span.worm=span)
	}

#	Enable REST-style online access to KEGG pathways
	if(!requireNamespace("KEGGREST",quietly=TRUE)) stop("KEGGREST required but is not available")

#	Convert to KEGG organism/species codes
	organism <- switch(species, "Hs"="hsa", "Dm"="dme", "Mm"="mmu", "Rn"="rno")

#	Get gene-KEGG mappings, and remove duplicate entries
	if(is.null(universe)) {

		path <- KEGGREST::keggLink("pathway", organism)
		path <- data.frame(kegg_id = names(path), path_id = path, stringsAsFactors = FALSE)
		if(organism == "dme") {
			EG <- KEGGREST::keggConv("dme", "ncbi-geneid")
			geneid <- names(EG)[match(path$kegg_id, EG)]
			geneid <- gsub("ncbi-geneid:", "", geneid)
		} else {
			geneid <- gsub(paste0(organism, ":"), "", path$kegg_id)
		}
		EG.KEGG <- data.frame(gene_id = geneid, path, stringsAsFactors = FALSE)

		universe <- unique(EG.KEGG$gene_id)
		universe <- as.character(universe)
	} else {

		universe <- as.character(universe)

		dup <- duplicated(universe)
		if(!is.null(prior.prob)) {
			if(length(prior.prob)!=length(universe)) stop("length(prior.prob) must equal length(universe)")
			prior.prob <- rowsum(prior.prob,group=universe,reorder=FALSE)
			n <- rowsum(rep_len(1L,length(universe)),group=universe,reorder=FALSE)
			prior.prob <- prior.prob/n
		}
		universe <- universe[!dup]

		path <- KEGGREST::keggLink("pathway", organism)
		path <- data.frame(kegg_id = names(path), path_id = path, stringsAsFactors = FALSE)
		if(organism == "dme") {
			EG <- KEGGREST::keggConv("dme", "ncbi-geneid")
			geneid <- names(EG)[match(path$kegg_id, EG)]
			geneid <- gsub("ncbi-geneid:", "", geneid)
		} else {
			geneid <- gsub(paste0(organism, ":"), "", path$kegg_id)
		}
		EG.KEGG <- data.frame(gene_id = geneid, path, stringsAsFactors = FALSE)

		m <- match(EG.KEGG$gene_id, universe)
		universe <- universe[m]
		if(!is.null(prior.prob)) prior.prob <- prior.prob[m]
		EG.KEGG <- EG.KEGG[EG.KEGG$gene_id %in% universe,]
	}

	Total <- length(unique(EG.KEGG$gene_id))
	if(Total<1L) stop("No genes found in universe")

#	Check prior probabilities
	if(!is.null(prior.prob)) {
		if(length(prior.prob)!=length(universe)) stop("length(prior.prob) must equal length(universe)")
	}

#	Overlap with DE genes
	isDE <- lapply(de, function(x) EG.KEGG$gene_id %in% x)
	TotalDE <- lapply(isDE, function(x) length(unique(EG.KEGG$gene_id[x])))
	nDE <- length(isDE)

	if(length(prior.prob)) {
	#	Probability weight for each gene
		m <- match(EG.KEGG$gene_id, universe)
		PW2 <- list(prior.prob[m])
		X <- do.call(cbind, c(N=1, isDE, PW=PW2))
	} else
		X <- do.call(cbind, c(N=1, isDE))

	group <- EG.KEGG$path_id
	S <- rowsum(X, group=group, reorder=FALSE)

	P <- matrix(0, nrow = nrow(S), ncol = nsets)

	if(length(prior.prob)) {

#		Calculate average prior prob for each set
		PW.ALL <- sum(prior.prob[universe %in% EG.KEGG$gene_id])
		AVE.PW <- S[,"PW"]/S[,"N"]
		W <- AVE.PW*(Total-S[,"N"])/(PW.ALL-S[,"N"]*AVE.PW)

#		Wallenius' noncentral hypergeometric test
		if(!requireNamespace("BiasedUrn",quietly=TRUE)) stop("BiasedUrn package required but is not available")
		for(j in 1:nsets) for(i in 1:nrow(S)) 
			P[i,j] <- BiasedUrn::pWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i],lower.tail=FALSE) + BiasedUrn::dWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i])
		S <- S[,-ncol(S)]

	} else {

#		Fisher's exact test
		for(j in 1:nsets)
			P[,j] <- phyper(q=S[,1+j]-0.5,m=TotalDE[[j]],n=Total-TotalDE[[j]], k=S[,"N"],lower.tail=FALSE)

	}

#	Assemble output
	g <- rownames(S)
	pathname <- KEGGREST::keggList("pathway")
	names(pathname) <- gsub("map", organism, names(pathname))
	m <- match(g, names(pathname))
	Results <- data.frame(Pathway = pathname[m], S, P, stringsAsFactors=FALSE)
	rownames(Results) <- g

#	Name p-value columns
	colnames(Results)[2+nsets+(1L:nsets)] <- paste0("P.", names(de))

	Results
}

topKEGG <- function(results, sort = NULL, number = 20L, truncate.path=NULL)
#	Extract top KEGG pathways from kegga output 
#	Gordon Smyth and Yifang Hu
#	Created 15 May 2015. Modified 4 August 2015.
{
#	Check results
	if(!is.data.frame(results)) stop("results should be a data.frame.")
	dimres <- dim(results)

#	Check number
	if(!is.numeric(number)) stop("number should be a positive integer")
	if(number > dimres[1L]) number <- dimres[1]
	if(number < 1L) return(results[integer(0),])

#	Number of gene lists for which results are reported
#	Lists are either called "Up" and "Down" or have user-supplied names
	nsets <- (dimres[2L]-2L) %/% 2L
	if(nsets < 1L) stop("results has wrong number of columns")
	setnames <- colnames(results)[3L:(2L+nsets)]

#	Check sort. Defaults to all gene lists.
	if(is.null(sort)) {
		isort <- 1L:nsets
	} else {
		sort <- as.character(sort)
		isort <- which(tolower(setnames) %in% tolower(sort))
		if(!length(isort)) stop("sort column not found in results")
	}

#	Sort by minimum p-value for specified gene lists
	P.col <- 2L+nsets+isort
	if(length(P.col)==1L) {
		o <- order(results[,P.col])
	} else {
		o <- order(do.call("pmin",as.data.frame(results[,P.col,drop=FALSE])))
	}
	tab <- results[o[1L:number],,drop=FALSE]

#	Truncate Pathway column for readability
	if(!is.null(truncate.path)) {
		truncate.path <- as.integer(truncate.path[1])
		truncate.path <- max(truncate.path,5L)
		truncate.path <- min(truncate.path,1000L)
		tm2 <- truncate.path-3L
		i <- (nchar(tab$Pathway) > tm2)
		tab$Pathway[i] <- paste0(substring(tab$Pathway[i],1L,tm2),"...")
	}

	tab
}
