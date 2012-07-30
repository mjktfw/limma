propTrueNull <- function(p)
#	Estimate proportion of null p-values
#	Belinda Phipson and Gordon Smyth
#	23 May 2012. Last revised 30 July 2012.
{
	n <- length(p)
	i <- n:1L
	p <- sort(p, decreasing = TRUE)
#	pi0 <- mean((1L:n)*p > 0.05)
	q <- pmin(n/i * p, 1)
	n1 <- n+1L
	sum(i*q) / n/n1*2
}
