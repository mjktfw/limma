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

setClass("exprSet2",representation(
	expressions="matrix",
	weights="matrix",
	targets="data.frame",
	probes="data.frame",
	printer="list",
	notes="character"
))

setAs("RGList", "exprSet2", function(from, to) {
	y <- new(to)
	if(length(from$G)) y@expressions <- cbind(from$G,from$R)
	if(length(from$weights)) y@weights <- cbind(from$weights,from$weights)
	if(length(from$genes)) y@probes <- from$genes
	if(length(from$printer)) y@printer <- unclass(from$printer)
	y
})

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
setIs("exprSet2","LargeDataObject")

setMethod("show","LargeDataObject",
#  Print and show method large data objects
#  Gordon Smyth
#  May 2003
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

dim.RGList <- function(x) dim(x$R)
dim.MAList <- function(x) dim(x$M)
dim.MArrayLM <- function(x) dim(x$coefficients)

as.MAList <- function(object) {
#	Convert marrayNorm object to MAList
#	Gordon Smyth
#	20 Sep 2003.  Last modified 20 Dec 2003.

	MA <- new("MAList")
	ifposlen <- function(x) if(length(x)) return(x) else return(NULL)
	MA$A <- ifposlen(object@maA)
	MA$M <- ifposlen(object@maM)
	MA$weights <- ifposlen(object@maW)
	MA$printer$ngrid.r <- ifposlen(object@maLayout@maNgr)
	MA$printer$ngrid.c <- ifposlen(object@maLayout@maNgc)
	MA$printer$nspot.r <- ifposlen(object@maLayout@maNsr)
	MA$printer$nspot.c <- ifposlen(object@maLayout@maNsc)
	MA$printer$notes <- ifposlen(object@maLayout@maNotes)
	MA$genes <- ifposlen(object@maGnames@maInfo)
	MA$genes$Labels <- ifposlen(object@maGnames@maLabels)
	attr(MA$genes,"notes") <- ifposlen(object@maGnames@maNotes)
	MA$genes$Sub <- ifposlen(object@maLayout@maSub)
	MA$genes$Plate <- ifposlen(object@maLayout@maPlate)
	MA$genes$Controls <- ifposlen(object@maLayout@maControls)
	MA$targets <- ifposlen(object@maTargets@maInfo)
	MA$targets$Labels <- ifposlen(object@maTargets@maLabels)
	MA$notes <- ifposlen(object@maNotes)
	MA$maNormCall <- ifposlen(object@maNormCall)
	MA
} 

