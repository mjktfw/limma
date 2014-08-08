goana <- function(de,...) UseMethod("goana")

goana.MArrayLM <- function(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, species = "Hs", trend = FALSE, ...)
#  Gene ontology analysis of DE genes from linear model fit
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014.  Last modified 23 July 2014.
{
	# Check fit
	if(is.null(de$p.value)) stop("p value not found in fit object (from eBayes).")	
	if(is.null(de$coefficients)) stop("coefficient not found in fit object.")
	if(length(coef) != 1) stop("coef length needs to be 1.")
	ngenes <- nrow(de)

	# Check geneid
	# Can be either a vector of IDs or a column name
	geneid <- as.character(geneid)
	if(length(geneid) == ngenes) {
		universe <- geneid
	} else
		if(length(geneid) == 1L) {
			universe <- de$genes[[geneid]]
			if(is.null(universe)) stop(paste("Column",geneid,"not found in de$genes"))
		} else
			stop("geneid has incorrect length")

	# Check trend
	# Can be logical, or a numeric vector of covariate values, or the name of the column containing the covariate values
	if(is.logical(trend)) {
		if(trend) {
			covariate <- de$Amean
			if(is.null(covariate)) stop("Amean not found in fit")
		}
	} else
		if(is.numeric(trend)) {
			if(length(trend) != ngenes) stop("If trend is numeric, then length must equal nrow(de)")
			covariate <- trend
			trend <- TRUE
		} else {
			if(is.character(trend)) {
				if(length(trend) != 1L) stop("If trend is character, then length must be 1")
				covariate <- de$genes[[trend]]
				if(is.null(covariate)) stop(paste("Column",trend,"not found in de$genes"))
				trend <- TRUE
			} else
				stop("trend is neither logical, numeric nor character")
		}

	# Check FDR
	if(!is.numeric(FDR) | length(FDR) != 1) stop("FDR must be numeric and of length 1.")
	if(FDR < 0 | FDR > 1) stop("FDR should be between 0 and 1.")

	# Get up and down DE genes
	fdr.coef <- p.adjust(de$p.value[,coef], method = "BH")
	EG.DE.UP <- universe[fdr.coef < FDR & de$coef[,coef] > 0]
	EG.DE.DN <- universe[fdr.coef < FDR & de$coef[,coef] < 0]
	de.gene <- list(Up=EG.DE.UP, Down=EG.DE.DN)

	# Fit monotonic cubic spline for DE genes vs. gene.weights
	if(trend) {
			PW <- isDE <- rep(0,nrow(de))
			isDE[fdr.coef < FDR] <- 1
			o <- order(covariate)
			PW[o] <- tricubeMovingAverage(isDE[o],span=0.5,full.length=TRUE)
	}
	if(!trend) PW <- NULL

	NextMethod(de = de.gene, universe = universe, species = species, weights = PW, ...)
}

goana.default <- function(de, universe = NULL, species = "Hs", weights = NULL, ...)
#  Gene ontology analysis of DE genes
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014.  Last modified 23 July 2014.
{
	# check de
	if(!is.list(de)) de <- list(DE1 = de)
	if(is.null(names(de))) names(de) <- paste0("DE", 1:length(de))
	de <- lapply(de, as.character)

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

	# Check universe
	if(is.null(universe)) universe <- EG.GO$gene_id
	universe <- as.character(universe)

	# Check weights
	if(!is.null(weights)){
		if(length(weights)!=length(universe)) stop("length(weights) must equal length(universe).")
	}

	# Reduce to universe
	EG.GO <- EG.GO[EG.GO$gene_id %in% universe, ]
	if(!length(EG.GO$gene_id)) stop("Universe is empty.")

	Total <- length(unique(EG.GO$gene_id))

	# Overlap with DE genes
	isDE <- lapply(de, function(x) EG.GO$gene_id %in% x)
	TotalDE <- lapply(isDE, function(x) length(unique(EG.GO$gene_id[x])))
	nDE <- length(isDE)

	if(length(weights)) {
		# Probability weight for each gene
		m <- match(EG.GO$gene_id, universe)
		PW2 <- list(weights[m])
		X <- do.call(cbind, c(N=1, isDE, PW=PW2))
	} else
		X <- do.call(cbind, c(N=1, isDE))

	group <- paste(EG.GO$go_id, EG.GO$Ontology, sep=".")
	S <- rowsum(X, group=group, reorder=FALSE)

	P <- matrix(0, nrow = nrow(S), ncol = nDE)

	if(length(weights)) {

		# Calculate weight
		require("BiasedUrn", character.only = TRUE)
		PW.ALL <- sum(weights[universe %in% EG.GO$gene_id])
		AVE.PW <- S[,"PW"]/S[,"N"]
		W <- AVE.PW*(Total-S[,"N"])/(PW.ALL-S[,"N"]*AVE.PW)

		# Wallenius' noncentral hypergeometric test
		for(j in 1:nDE){

			for(i in 1:nrow(S)){

				P[i,j] <- pWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i],lower.tail=FALSE)+
					dWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i])
			}

		}

		S <- S[,-ncol(S)]

	} else {

		# Fisher's exact test
		for(j in 1:nDE){

			P[,j] <- phyper(q=S[,1+j]-0.5,m=TotalDE[[j]],n=Total-TotalDE[[j]], k=S[,"N"],lower.tail=FALSE)
		}

	}

	# Assemble output
	g <- strsplit2(rownames(S),split="\\.")
	TERM <- select(GO.db,keys=g[,1],columns="TERM")
	Results <- data.frame(Term = TERM[[2]], Ont = g[,2], S, P, stringsAsFactors=FALSE)
	rownames(Results) <- g[,1]

	# Name P value for the DE genes
	iTON <- c(1:3)
	iDE <- 3+c(1:nDE)
	PDE<- paste0("P.", colnames(Results)[iDE])
	colnames(Results)[-c(iTON,iDE)] <- PDE

	Results
}

topGO <- function(results, ontology = c("BP", "CC", "MF"), sort = "up", number = 20L){
#  Extract sorted goana gene ontology test results 
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014. Last modified 21 July 2014.

	# Check results
	if(!is.data.frame(results)) stop("Expect a dataframe with goana results.")

	# Check ontology
	ontology <- match.arg(ontology, c("BP", "CC", "MF"), several.ok = TRUE)

	# Check sort and sort by P value
	if(length(sort) != 1) stop("sort length needs to be 1.")
	sort <- tolower(colnames(results)) %in% tolower(paste0("P.", sort))
	if(sum(sort)!=1) stop("sort not found.")

	# Check number
	if(!is.numeric(number)) stop("Need to input number.")
	if(number < 1L) return(data.frame())

	# Select ontologies
	sel <- results$Ont %in% ontology
	results <- results[sel,]

	# Sort by p value
	o <- order(results[, sort], rownames(results))
	results <- results[o,]

	if(number >= nrow(results)) results else results[1:number,]
}