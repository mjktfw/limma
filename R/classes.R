#	CLASSES.R

setClass("RGList",
#  Class to hold initial read-in data
representation("list")
)

setClass("MAList",
#  Class to hold normalized, rotated data
representation("list")
)

setClass("MArrayLM",
#  Linear model fit
representation("list")
)

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

