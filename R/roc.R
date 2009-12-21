#  ROC.R

auROC <- function(truth, stat=NULL) {
#	Area under Receiver Operating Curve for empirical data
#	Gordon Smyth
#	21 Dec 2003. Last modified 4 Dec 2009.

	if(!length(truth)) return(NULL)
	truth <- as.numeric(as.logical(truth))
	if(!is.null(stat)) {
		if(length(stat) != length(truth)) stop("lengths differ")
		truth[is.na(stat)] <- NA
		truth <- truth[order(stat,decreasing=TRUE)]
	}
	isna <- is.na(truth)
	if(any(isna)) truth <- truth[!isna]
	sens <- cumsum(truth)/sum(truth)
	mean(sens[truth==0])
}

