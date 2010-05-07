##  ALIAS2SYMBOL.R

alias2Symbol <- function(alias,species="Hs",expand.symbols=FALSE)
#  Convert a set of alias names to official gene symbols
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  4 Sep 2008. Last revised 15 Jan 2009
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))
	DB <- paste("org",species,"eg","db",sep=".")
	ALIAS2EG <- paste("org",species,"egALIAS2EG",sep=".")
	SYMBOL <- paste("org",species,"egSYMBOL",sep=".")
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
#  Convert a set of alias names to official gene symbols of the same length
#  via Entrez Gene identifiers
#  Di Wu, Gordon Smyth and Yifang Hu
#  3 Sep 2009.  Last modified 17 Dec 2009.
{
	alias <- as.character(alias)
	species <- match.arg(species,c("Dm","Hs","Mm","Rn"))
	DB <- paste("org",species,"eg","db",sep=".")
	ALIAS2EG <- paste("org",species,"egALIAS2EG",sep=".")
	SYMBOL <- paste("org",species,"egSYMBOL",sep=".")
	suppressPackageStartupMessages(require(DB,character.only=TRUE))
	
	isSymbol <- alias %in% Rkeys(get(SYMBOL)) 
	Symbol<-rep.int(NA,length(alias))
	Symbol[isSymbol]<-alias[isSymbol]
				
	isalias<-(alias[!isSymbol]) %in% (Rkeys(get(ALIAS2EG)))
	alias2<-(alias[!isSymbol])[isalias]

	aliasTbl<-toTable(get(ALIAS2EG)[alias2])
	hits<-names(table(aliasTbl$alias_symbol))[as.numeric(table(aliasTbl$alias_symbol))>1]
	if(length(hits)>0) warning("Multiple Hits for ", hits)
	
	aliasTbl.o<-aliasTbl[match(alias2,aliasTbl$alias_symbol),]
	symb<-toTable(get(SYMBOL)[aliasTbl.o$gene_id])
	Symbol[!isSymbol][isalias]<-symb[match(aliasTbl.o$gene_id,symb$gene_id),]$symbol
	Symbol	
}

