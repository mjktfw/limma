#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG) {
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003. Last modified 26 April 2003.

	RG$R <- RG$R - RG$Rb
	RG$G <- RG$G - RG$Gb
	RG$Rb <- NULL
	RG$Gb <- NULL
	new("RGList",unclass(RG))
}
