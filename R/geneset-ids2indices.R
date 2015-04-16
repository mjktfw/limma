ids2indices <- function(gene.sets, identifiers, remove.empty=TRUE)
# Make a list of gene identifiers into a list of indices for gene sets
# Gordon Smyth and Yifang Hu
# 25 March 2009.  Last modified 19 June 2014.
{
	gene.sets <- as.list(gene.sets)
	index <- lapply(gene.sets, function(x) which(identifiers %in% x))
	if(remove.empty)
		for (i in length(index):1) if(!length(index[[i]])) index[[i]] <- NULL
	index
}
