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
		isDE <- as.numeric(fdr.coef < FDR)
		o <- order(covariate)
		PW <- rep(0,nrow(de))
		PW[o] <- tricubeMovingAverage(isDE[o],span=0.5,full.length=TRUE)
	}
	if(!trend) PW <- NULL

	NextMethod(de = de.gene, universe = universe, species = species, prior.prob = PW, ...)
}

goana.default <- function(de, universe = NULL, species = "Hs", prior.prob = NULL, ...)
#  Gene ontology analysis of DE genes
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014.  Last modified 28 August 2014.
{
	# Ensure de is a list
	if(!is.list(de)) de <- list(DE = de)

	# Stop if components of de are not vectors
	if(!all(vapply(de,is.vector,TRUE))) stop("components of de should be vectors")

	# Ensure gene IDs are of character mode
	de <- lapply(de, as.character)

	# Ensure all gene sets have names
	nsets <- length(de)

	names(de) <- trimWhiteSpace(names(de))
	NAME <- names(de)
	i <- NAME == "" | is.na(NAME)
	NAME[i] <- paste0("DE", 1:nsets)[i]
	names(de) <- makeUnique(NAME)

	# Select species
	species <- match.arg(species, c("Hs", "Mm", "Rn", "Dm"))

	# Load package of GO terms
	if(!suppressPackageStartupMessages(require("GO.db", character.only = TRUE))) stop("Go.db package not available")

	# Load species annotation package
	DB <- paste("org", species, "eg", "db", sep = ".")
	if(!suppressPackageStartupMessages(require(DB, character.only = TRUE))) stop(DB,"package not available")

	# Get gene-GOterm mappings, and remove duplicate entries
	GO2ALLEGS <- paste("org", species, "egGO2ALLEGS", sep = ".")
	if(is.null(universe)) {
		EG.GO <- toTable(get(GO2ALLEGS))
		d <- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
		EG.GO <- EG.GO[!d, ]
		universe <- unique(EG.GO$gene_id)
		universe <- as.character(universe)
	} else {

		universe <- as.character(universe)

		dup <- duplicated(universe)
		if(!is.null(prior.prob)) {
			if(length(prior.prob)!=length(universe)) stop("length(prior.prob) must equal length(universe)")
			prior.prob <- prior.prob[!dup]
		}
		universe <- universe[!dup]

		GO2ALLEGS <- get(GO2ALLEGS)
		m <- match(Lkeys(GO2ALLEGS),universe,0L)
		universe <- universe[m]
		if(!is.null(prior.prob)) prior.prob <- prior.prob[m]

		Lkeys(GO2ALLEGS) <- universe
		EG.GO <- toTable(GO2ALLEGS)
		d <- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
		EG.GO <- EG.GO[!d, ]
	}

	Total <- length(unique(EG.GO$gene_id))
	if(Total<1L) stop("No genes found in universe")

	# Check prior probabilities
	if(!is.null(prior.prob)) {
		if(length(prior.prob)!=length(universe)) stop("length(prior.prob) must equal length(universe)")
	}

	# Overlap with DE genes
	isDE <- lapply(de, function(x) EG.GO$gene_id %in% x)
	TotalDE <- lapply(isDE, function(x) length(unique(EG.GO$gene_id[x])))
	nDE <- length(isDE)

	if(length(prior.prob)) {
		# Probability weight for each gene
		m <- match(EG.GO$gene_id, universe)
		PW2 <- list(prior.prob[m])
		X <- do.call(cbind, c(N=1, isDE, PW=PW2))
	} else
		X <- do.call(cbind, c(N=1, isDE))

	group <- paste(EG.GO$go_id, EG.GO$Ontology, sep=".")
	S <- rowsum(X, group=group, reorder=FALSE)

	P <- matrix(0, nrow = nrow(S), ncol = nsets)

	if(length(prior.prob)) {

		# Calculate weight
		require("BiasedUrn", character.only = TRUE)
		PW.ALL <- sum(prior.prob[universe %in% EG.GO$gene_id])
		AVE.PW <- S[,"PW"]/S[,"N"]
		W <- AVE.PW*(Total-S[,"N"])/(PW.ALL-S[,"N"]*AVE.PW)

		# Wallenius' noncentral hypergeometric test
		for(j in 1:nsets) for(i in 1:nrow(S))
			P[i,j] <- pWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i],lower.tail=FALSE) + dWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i])

		S <- S[,-ncol(S)]

	} else {

		# Fisher's exact test
		for(j in 1:nsets)
			P[,j] <- phyper(q=S[,1+j]-0.5,m=TotalDE[[j]],n=Total-TotalDE[[j]], k=S[,"N"],lower.tail=FALSE)

	}

	# Assemble output
	g <- strsplit2(rownames(S),split="\\.")
	TERM <- select(GO.db,keys=g[,1],columns="TERM")
	Results <- data.frame(Term = TERM[[2]], Ont = g[,2], S, P, stringsAsFactors=FALSE)
	rownames(Results) <- g[,1]

	# Name p-value columns
	colnames(Results)[3+nsets+(1L:nsets)] <- paste0("P.", names(de))

	Results
}

topGO <- function(results, ontology = c("BP", "CC", "MF"), sort = NULL, number = 20L)
#  Extract sorted results from goana output 
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014. Last modified 29 August 2014.
{
	# Check results
	if(!is.data.frame(results)) stop("results should be a data.frame.")

	# Number of gene sets
	nsets <- (ncol(results)-3L) %/% 2L
	setnames <- colnames(results)[4L:(3L+nsets)]

	# Check ontology
	ontology <- match.arg(unique(ontology), c("BP", "CC", "MF"), several.ok = TRUE)

	# Check sort and find p-value column
	if(is.null(sort)) {
		P.col <- 4L+nsets
	} else {
		sort <- as.character(sort[1])
		if(sort=="up") sort="Up"
		if(sort=="down") sort="Down"
		sort <- paste0("^",sort,"$")
		P.col <- grep(sort, setnames)
		if(!length(P.col)) stop("sort column not found")
		P.col <- 3L+nsets+P.col
	}

	# Check number
	if(!is.numeric(number)) stop("Need to input number.")
	if(number < 1L) return(data.frame())

	# Limit results to specified ontologies
	if(length(ontology) < 3L) {
		sel <- results$Ont %in% ontology
		results <- results[sel,]
	}

	# Sort by p-value
	o <- order(results[,P.col], rownames(results))
	if(number < length(o)) o <- o[1:number]

	results[o,,drop=FALSE]
}
