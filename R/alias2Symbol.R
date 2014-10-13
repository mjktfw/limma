##  ALIAS2SYMBOL.R

alias2Symbol <- function(alias,species="Hs",expand.symbols=FALSE)
#  Convert a set of alias names to official gene symbols
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  4 Sep 2008. Last revised 29 September 2013
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))
	DB <- paste("org",species,"eg","db",sep=".")
	ALIAS2EG <- paste("org",species,"egALIAS2EG",sep=".")
	SYMBOL <- paste("org",species,"egSYMBOL",sep=".")
	suppressPackageStartupMessages(require("AnnotationDbi",character.only=TRUE))
	suppressPackageStartupMessages(require(DB,character.only=TRUE))
	if(expand.symbols)
	{
		alias <- intersect(alias,Rkeys(get(ALIAS2EG)))
		eg <- mappedLkeys(get(ALIAS2EG)[alias])
		mappedRkeys(get(SYMBOL)[eg])
	}
	else
	{
		isSymbol <- alias %in% Rkeys(get(SYMBOL)) 
		alias2 <- intersect(alias[!isSymbol],Rkeys(get(ALIAS2EG)))
		eg <- mappedLkeys(get(ALIAS2EG)[alias2])
		c(alias[isSymbol],mappedRkeys(get(SYMBOL)[eg]))

	}
}

alias2SymbolTable <- function(alias,species="Hs")
#  Convert a vector of alias names to the vector of corresponding official gene symbols
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  Created 3 Sep 2009.  Last modified 29 September Dec 2013.
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))
	DB <- paste("org",species,"eg","db",sep=".")
	ALIAS2EG <- paste("org",species,"egALIAS2EG",sep=".")
	SYMBOL <- paste("org",species,"egSYMBOL",sep=".")
	suppressPackageStartupMessages(require("AnnotationDbi",character.only=TRUE))
	suppressPackageStartupMessages(require(DB,character.only=TRUE))

	isSymbol <- alias %in% Rkeys(get(SYMBOL))
	Symbol <- alias
	Symbol[!isSymbol] <- NA

	OtherAliases <- alias[!isSymbol]
	isAlias <- OtherAliases %in% Rkeys(get(ALIAS2EG))
	if(!any(isAlias)) return(Symbol)
	OtherAliases <- OtherAliases[isAlias]

	AliasTbl <- toTable(get(ALIAS2EG)[OtherAliases])
	if(anyDuplicated(AliasTbl$alias_symbol)) warning("Multiple symbols ignored for one or more aliases")
	SymbolTbl <- toTable(get(SYMBOL)[AliasTbl$gene_id])
	m <- match(OtherAliases,AliasTbl$alias_symbol)
	GeneID <- AliasTbl$gene_id[m]
	m <- match(GeneID,SymbolTbl$gene_id)
	Symbol[!isSymbol][isAlias] <- SymbolTbl$symbol[m]
	Symbol
}

