tricubeMovingAverage <- function(x,span=0.5,full.length=TRUE) {
	n <- length(x)
	width <- span*n
#	Round width of smoothing window to nearest odd number
	hwidth <- as.integer(width %/% 2L)
	if(hwidth <= 0L) return(x)
	width <- 2L * hwidth + 1L
	u <- seq(from=-1,to=1,length=width) * width / (width+1)
	tricube.weights <- (1-abs(u)^3)^3
	tricube.weights <- tricube.weights/sum(tricube.weights)
	if(!full.length) return(as.vector(filter(x,tricube.weights),mode="numeric")[(hwidth+1):(n-hwidth)])
	z <- numeric(hwidth)
	x <- as.vector(filter(c(z,x,z),tricube.weights),mode="numeric")[(hwidth+1):(n+hwidth)]
	cw <- cumsum(tricube.weights)
	x[1:hwidth] <- x[1:hwidth] / cw[(width-hwidth):(width-1)]
	x[(n-hwidth+1):n] <- x[(n-hwidth+1):n] / cw[(width-1):(width-hwidth)]
	x
}
