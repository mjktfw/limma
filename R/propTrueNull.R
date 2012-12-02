propTrueNull <- function(p, method="lfdr", nbins=20, ...)
#	Estimate proportion of null p-values
#	Belinda Phipson and Gordon Smyth
#	Created 23 May 2012. Last revised 2 Dec 2012.
{
	method <- match.arg(method, c("lfdr","hist","convest"))
	switch(method,
		lfdr = .propTrueNullByLocalFDR(p),
		hist = .propTrueNullFromHistogram(p, nbins=nbins),
		convest = convest(p, ...)
	)
}

.propTrueNullByLocalFDR <- function(p)
#	Estimate proportion of null p-values
#	by average local FDR
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

.propTrueNullFromHistogram <- function(p, nbins=20) 
#	Estimate proportion of null p-values
#	by histogram method

#	Adapted by Gordon Smyth from the function estimate.m0 from
#	http://www.public.iastate.edu/~dnett/microarray/multtest.txt
#	Accessed March 2012

#	Created 2 Dec 2012. Last revised 2 Dec 2012.
{
	bin <- c(-0.1, (1:nbins)/nbins)
	bin.counts <- tabulate(cut(p,bin))
	tail.means <- rev(cumsum(rev(bin.counts))/(1:nbins))
	index <- which(tail.means >= bin.counts)[1]
	tail.means[index]/tail.means[1]
}


