#	CLASSES.R

setClass("RGList",
#  Class to hold initial read-in data
representation("list")
)

#setValidity("RGList",
#function(object) {
#	R <- object$R
#	G <- object$G
#	if(is.null(R) || is.null(G)) return("Element R or G missing")
#	if(!identical(dim(R),dim(G))) return("Dimensions of R and G don't match")
#	if(!isNumeric(R) || !isNumeric(G)) return("R or G contain non-numeric elements")
#	if(length(dim(R)) > 2) return("R and G have more than two dimensions")
#})

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

#  Linear model fit
setClass("MArrayLM", representation(
	genes="data.frame",
	design="matrix",
	contrasts="matrix",
	coefficients="matrix",
	stdev.unscaled="matrix",
	s2.residual="numeric",
	df.residual="numeric",
	s2.prior="numeric",
	df.prior="numeric",
	s2.post="numeric",
	tstat="matrix",
	varcoef.prior="numeric",
	call="call"
))

printHead <- function(x) {
	what <- "other"
	if(is.vector(x)) what <- "vector"
	if(is.matrix(x) || is.data.frame(x)) what <- "TwoD"
	switch(what,
		vector={
			n <- length(x)
			if(n > 20) {
				print(x[1:5])
				cat(n-5,"more elements ...\n")
			} else
				print(x)
		},
		TwoD={
			n <- nrow(x)
			if(n > 10) {
				print(x[1:5,])
				cat(n-5,"more rows ...\n")
			} else
				print(x)
		},
		other=print(x)
	)
}

setClass("LargeDataObject")
setIs("RGList","LargeDataObject")
setIs("MAList","LargeDataObject")
setIs("MArrayLM","LargeDataObject")

setMethod("show","LargeDataObject",
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
	for (what in names(object)) {
		x <- object[[what]]
		cat("$",what,"\n",sep="")
		printHead(x)
		cat("\n")
	}
	for (what in setdiff(slotNames(object),".Data")) {
		x <- slot(object,what)
		if(length(x) > 0) {
			cat("@",what,"\n",sep="")
			printHead(x)
			cat("\n")
		}
	}
})
