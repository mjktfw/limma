#  UTILITY FUNCTIONS

matvec <- function(M,v) {
#	Multiply the columns of matrix by the elements of a vector,
#	i.e., compute M %*% diag(v)
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[2]) stop("Dimensions do not match")
	t(v * t(M))
}

vecmat <- function(v,M) {
#	Multiply the rows of matrix by the elements of a vector,
#	i.e., compute diag(v) %*% M
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[1]) stop("Dimensions do not match")
	v * M
}

isNumeric <- function(x) {
#	Test for numeric argument or data.frame with numeric columns
#	Gordon Smyth
#	12 April 2003

	is.numeric(x) || (is.data.frame(x) && length(x)>0 && all(unlist(lapply(x,is.numeric))))
}

helpMethods <- function(genericFunction) {
#	Prompt user for help topics on methods for generic function
#	Gordon Smyth
#	21 April 2003

	objectclass <- class(genericFunction)
 	if(objectclass != "genericFunction") {
		if(objectclass == "character" && isGeneric(genericFunction))
			genericFunction <- getGeneric(genericFunction)
		else {
			cat("Not a generic function\n")
			return(invisible())
		}
	}
	functionname <- genericFunction@generic
	methodnames <- names(getMethods(genericFunction)@methods)
	nmethods <- length(methodnames)
	if(nmethods == 0) {
		cat("No available methods\n")
		return(invisible())
	}
	aliasnames <- paste(functionname,methodnames,sep=".")
	for (i in 1:nmethods) cat(i,": ",aliasnames[i],"\n",sep="")
	cat("Type number to choose help topic: ")
	n <- as.integer(readline())
	if(n > 0 && n <= nmethods)
		eval(parse(text=paste("help(",aliasnames[n],")",sep="")))
	else {
	 	cat("No topic chosen\n")
	 	return(invisible())
	}
}
