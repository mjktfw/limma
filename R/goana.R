goana <- function(fit, coef = ncol(fit), geneid = "GeneID", FDR = 0.05, species = "Hs"){
#  Gene ontology analysis of DE genes from linear model fit
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014. Last modified 24 June 2014.

	# Check input
	if(!is(fit, "MArrayLM")) stop("fit must be an MArrayLM object.")
	if(is.null(fit$p.value)) stop("p value not found in fit object (from eBayes).")
	if(is.null(fit$coefficients)) stop("coefficient not found in fit object.")

	if(length(geneid) == nrow(fit)){

		fit$genes$GeneID <- geneid
		EG.col <- "GeneID"
	}

	else if(length(geneid) == 1) EG.col <- as.character(geneid)

	if(is.null(fit$genes)) stop("no annotation (genes) found.")

	EG.All <- as.character(fit$genes[[EG.col]])

	# Get up and down DE genes
	fdr.coef <- p.adjust(fit$p.value[,coef], method = "BH")

	EG.DE.UP <- EG.All[fdr.coef < FDR & fit$coef[,coef] > 0]
	EG.DE.DN <- EG.All[fdr.coef < FDR & fit$coef[,coef] < 0]

	# Select species
	species <- match.arg(species, c("Hs", "Mm", "Rn", "Dm"))

	# Load packages
	DB <- paste("org", species, "eg", "db", sep = ".")
	suppressPackageStartupMessages(require(DB, character.only = TRUE))
	suppressPackageStartupMessages(require("GO.db", character.only = TRUE))

	# Get GO data
	GO2ALLEGS <- paste("org", species, "egGO2ALLEGS", sep = ".")
	EG.GO <- toTable(get(GO2ALLEGS))

	# Remove duplicates
	d <- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
	EG.GO <- EG.GO[!d, ]

	# Reduce to universe
	EG.GO <- EG.GO[EG.GO$gene_id %in% EG.All, ]
	Total <- length(unique(EG.GO$gene_id))

	# Overlap with DE genes
	isDE.UP <- (EG.GO$gene_id %in% EG.DE.UP)
	isDE.DN <- (EG.GO$gene_id %in% EG.DE.DN)
	TotalDE.UP <- length(unique(EG.GO$gene_id[isDE.UP]))
	TotalDE.DN <- length(unique(EG.GO$gene_id[isDE.DN]))

	X <- cbind(N=1, Up=isDE.UP, Down=isDE.DN)
	group <- paste(EG.GO$go_id, EG.GO$Ontology, sep=".")
	S <- rowsum(X, group=group, reorder=FALSE)

	# Fisher's exact test
	p.UP <- phyper(q=S[,"Up"]-0.5,m=TotalDE.UP,n=Total-TotalDE.UP,k=S[,"N"],lower.tail=FALSE)
	p.DN <- phyper(q=S[,"Down"]-0.5,m=TotalDE.DN,n=Total-TotalDE.DN,k=S[,"N"],lower.tail=FALSE)

	# Assemble output
	g <- strsplit2(rownames(S),split="\\.")
	TERM <- select(GO.db,keys=g[,1],columns="TERM")
	Results <- data.frame(Term=TERM[[2]],Ont=g[,2],S,P.Up=p.UP,P.Down=p.DN,stringsAsFactors=FALSE)
	rownames(Results) <- g[,1]

	Results
}

topGO <- function(results, ontology = c("BP", "CC", "MF"), sort = "up", number = 20L){
#  Extract sorted goana gene ontology test results 
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014. Last modified 24 June 2014.
	
	# Check input
	if(any(! c("Term", "N", "Up", "Down", "Ont", "P.Up", "P.Down") %in% colnames(results))) stop("Results column names don't match GOTest results.")
	
	ontology <- match.arg(ontology, c("BP", "CC", "MF"), several.ok = TRUE)
	
	sort <- match.arg(tolower(sort), c("up", "down"))
	sort <- switch(sort, up = "P.Up", down = "P.Down")
	
	if(!is.numeric(number)) stop("Need to input number.")
	if(number < 1L) return(data.frame())

	# Select ontologies
	sel <- results$Ont %in% ontology
	results <- results[sel,]

	# Sort by Up or Down p value
	o <- order(results[, sort], rownames(results))
	results <- results[o,]

	if(number >= nrow(results)) results else results[1:number,]
}

