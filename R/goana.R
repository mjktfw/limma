goana <- function(de,...) UseMethod("goana")

goana.MArrayLM <- function(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, species = "Hs", trend = FALSE, plot=FALSE, ...)
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

	# Fit trend in DE with respect to the covariate
	if(trend) {
		isDE <- as.numeric(fdr.coef < FDR)
		o <- order(covariate)
		PW <- rep(0,nrow(de))
		PW[o] <- tricubeMovingAverage(isDE[o],span=0.5,full.length=TRUE)
		if(plot) plot(covariate[o],PW[o],type="l",xlab="Abundance",ylab="prior.prob")
	}
	if(!trend) PW <- NULL

	goana(de = de.gene, universe = universe, species = species, prior.prob = PW, ...)
}

goana.default <- function(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL, plot=FALSE, ...)
#  Gene ontology analysis of DE genes
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014.  Last modified 27 March 2015.
{
	# Ensure de is a list
	if(!is.list(de)) de <- list(DE = de)

	# Stop if components of de are not vectors
	if(!all(vapply(de,is.vector,TRUE))) stop("components of de should be vectors")

	# Ensure gene IDs are of character mode
	de <- lapply(de, as.character)
	if(!is.null(universe)) universe <- as.character(universe)

	# Ensure all gene sets have unique names
	nsets <- length(de)
	names(de) <- trimWhiteSpace(names(de))
	NAME <- names(de)
	i <- which(NAME == "" | is.na(NAME))
	NAME[i] <- paste0("DE",i)
	names(de) <- makeUnique(NAME)

	# Select species
	species <- match.arg(species, c("Hs", "Mm", "Rn", "Dm"))

	# Fit trend in DE with respect to the covariate, combining all de lists
	if(!is.null(covariate)) {
		isDE <- as.numeric(universe %in% unlist(de))
		o <- order(covariate)
		prior.prob <- covariate
		prior.prob[o] <- tricubeMovingAverage(isDE[o],span=0.5,full.length=TRUE)
		if(plot) plot(covariate[o],prior.prob[o],type="l")
	}

	# Load package of GO terms
	if(!suppressPackageStartupMessages(require("GO.db", character.only = TRUE))) stop("GO.db package not available")

	# Load species annotation package
	DB <- paste("org", species, "eg", "db", sep = ".")
	if(!suppressPackageStartupMessages(require(DB, character.only = TRUE))) stop(DB,"package not available")

	# Get gene-GOterm mappings, and remove duplicate entries
	GO2ALLEGS <- paste("org", species, "egGO2ALLEGS", sep = ".")
	if(is.null(universe)) {
		EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
		d <- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
		EG.GO <- EG.GO[!d, ]
		universe <- unique(EG.GO$gene_id)
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

		GO2ALLEGS <- get(GO2ALLEGS)
		m <- match(AnnotationDbi::Lkeys(GO2ALLEGS),universe,0L)
		universe <- universe[m]
		if(!is.null(prior.prob)) prior.prob <- prior.prob[m]

		AnnotationDbi::Lkeys(GO2ALLEGS) <- universe
		EG.GO <- AnnotationDbi::toTable(GO2ALLEGS)
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
		if(!requireNamespace("BiasedUrn",quietly=TRUE)) stop("BiasedUrn package required but is not available")
		PW.ALL <- sum(prior.prob[universe %in% EG.GO$gene_id])
		AVE.PW <- S[,"PW"]/S[,"N"]
		W <- AVE.PW*(Total-S[,"N"])/(PW.ALL-S[,"N"]*AVE.PW)

		# Wallenius' noncentral hypergeometric test
		for(j in 1:nsets) for(i in 1:nrow(S))
			P[i,j] <- BiasedUrn::pWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i],lower.tail=FALSE) + BiasedUrn::dWNCHypergeo(S[i,1+j], S[i,"N"], Total-S[i,"N"], TotalDE[[j]], W[i])

		S <- S[,-ncol(S)]

	} else {

		# Fisher's exact test
		for(j in 1:nsets)
			P[,j] <- phyper(q=S[,1+j]-0.5,m=TotalDE[[j]],n=Total-TotalDE[[j]], k=S[,"N"],lower.tail=FALSE)

	}

	# Assemble output
	g <- strsplit2(rownames(S),split="\\.")
	TERM <- AnnotationDbi::select(GO.db::GO.db,keys=g[,1],columns="TERM")
	Results <- data.frame(Term = TERM[[2]], Ont = g[,2], S, P, stringsAsFactors=FALSE)
	rownames(Results) <- g[,1]

	# Name p-value columns
	colnames(Results)[3+nsets+(1L:nsets)] <- paste0("P.", names(de))

	Results
}

topGO <- function(results, ontology = c("BP", "CC", "MF"), sort = NULL, number = 20L)
#  Extract top GO terms from goana output 
#  Gordon Smyth and Yifang Hu
#  Created 20 June 2014. Last modified 4 April 2015.
{
	# Check results
	if(!is.data.frame(results)) stop("results should be a data.frame.")

	# Check ontology
	ontology <- match.arg(unique(ontology), c("BP", "CC", "MF"), several.ok = TRUE)

	# Limit results to specified ontologies
	if(length(ontology) < 3L) {
		sel <- results$Ont %in% ontology
		results <- results[sel,]
	}
	dimres <- dim(results)

	# Check number
	if(!is.numeric(number)) stop("number should be a positive integer")
	if(number > dimres[1L]) number <- dimres[1]
	if(number < 1L) return(results[integer(0),])

	# Number of gene lists for which results are reported
	# Lists are either called "Up" and "Down" or have user-supplied names
	nsets <- (dimres[2L]-3L) %/% 2L
	if(nsets < 1L) stop("results has wrong number of columns")
	setnames <- colnames(results)[4L:(3L+nsets)]

	# Check sort. Defaults to all gene lists.
	if(is.null(sort)) {
		isort <- 1L:nsets
	} else {
		sort <- as.character(sort)
		isort <- which(tolower(setnames) %in% tolower(sort))
		if(!length(isort)) stop("sort column not found in results")
	}

	# Sort by minimum p-value for specified gene lists
	P.col <- 3L+nsets+isort
	if(length(P.col)==1L) {
		o <- order(results[,P.col])
	} else {
		o <- order(do.call("pmin",as.data.frame(results[,P.col,drop=FALSE])))
	}
	results[o[1L:number],,drop=FALSE]
}
