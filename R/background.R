#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG, method="subtract") {
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003.  Last modified 22 August 2003.

	method <- match.arg(method, c("subtract", "positive"))
	RG$R <- RG$R - RG$Rb
	RG$G <- RG$G - RG$Gb
	if(method=="positive") {
		RG$R[ RG$R==0 ] <- 0.75
		RG$R[ RG$G==0 ] <- 0.75
		RG$R[ RG$R<0 ] <- 0.5
		RG$G[ RG$G<0 ] <- 0.5
	}
	RG$Rb <- NULL
	RG$Gb <- NULL
	new("RGList",unclass(RG))
}
