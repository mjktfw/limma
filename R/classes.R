#	CLASSES.R

setClass("RGList",
#  Class to hold initial read-in data
representation("list")
)

setValidity("RGList",
function(object) {
	R <- object$R
	G <- object$G
	if(is.null(R) || is.null(G)) return("Element R or G missing")
	if(!identical(dim(R),dim(G))) return("Dimensions of R and G don't match")
	if(!isNumeric(R) || !isNumeric(G)) return("R or G contain non-numeric elements")
	if(length(dim(R)) > 2) return("R and G have more than two dimensions")
})

setClass("MAList",
#  Class to hold normalized, rotated data
representation("list")
)

##  Matrices of log expression values
#setClass("ExpressionMatrix", representation(
#	"matrix",
#	background="matrix",
#	weights="matrix",
#	printlayout="list",
#	replicatespots="list"
#))

##  Linear model fit
#setClass("MArrayLM", representation(
#	coefficients="matrix",
#	stdev.unscaled="matrix",
#	s2.residual="numeric",
#	df.residual="numeric",
#	correlation="numeric",
#	s2.prior="numeric",
#	df.prior="numeric",
#	var.prior="numeric"
#))
