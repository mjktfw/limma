#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG, method="subtract") {
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003.  Last modified 6 September 2003.

	method <- match.arg(method, c("subtract", "half", "minimum"))
	RG$R <- RG$R - RG$Rb
	RG$G <- RG$G - RG$Gb
	if(method=="half") {
		RG$R <- pmax(RG$R, 0.5)
		RG$G <- pmax(RG$G, 0.5)
	}
	if(method=="minimum") {
		RG$R <- as.matrix(RG$R)
		RG$G <- as.matrix(RG$G)
		for (slide in 1:ncol(RG$R)) {
			i <- RG$R[,slide] < 1e-18
			if(any(i)) {
				m <- min(RG$R[!i,slide],na.rm=TRUE)
				RG$R[i,slide] <- m/2
			}
			i <- RG$G[,slide] < 1e-18
			if(any(i)) {
				m <- min(RG$G[!i,slide],na.rm=TRUE)
				RG$G[i,slide] <- m/2
			}
		}
	}
	RG$Rb <- NULL
	RG$Gb <- NULL
	new("RGList",unclass(RG))
}

