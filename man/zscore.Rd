\name{zscore}
\alias{zscore}
\alias{zscoreGamma}
\alias{zscoreT}
\alias{tZscore}
\alias{zscoreHyper}

\title{Z-score Equivalents}

\description{
Compute z-score equivalents of non-normal random deviates.
}

\usage{
zscore(q, distribution, ...)
zscoreGamma(q, shape, rate = 1, scale = 1/rate) 
zscoreT(x, df)
tZscore(x, df)
zscoreHyper(q, m, n, k) 
}

\arguments{
\item{q, x}{numeric vector or matrix giving deviates of a random variable}
\item{distribution}{character name of probabability distribution for which a cumulative distribution function exists}
\item{\ldots}{other arguments specify distributional parameters and are passed to the cumulative distribution}
\item{shape}{gamma shape parameter (>0)}
\item{rate}{gamma rate parameter (>0)}
\item{scale}{gamma scale parameter (>0)}
\item{df}{degrees of freedom (>0 for \code{zscoreT} or >=1 for \code{tZscore})}
\item{m}{as for \code{\link{qhyper}}}
\item{n}{as for \code{\link{qhyper}}}
\item{k}{as for \code{\link{qhyper}}}
}

\value{
Numeric vector giving equivalent deviates from the standard normal distribution.
The exception is \code{tZscore} which gives deviates from the specified t-distribution.
}

\details{
These functions compute the standard normal deviates which have the same quantiles as the given values in the specified distribution.
For example, if \code{z <- zscoreT(x,df=df)} then \code{pnorm(z)} equals \code{pt(x,df=df)}.

\code{zscore} works for any distribution for which a cumulative distribution function (like \code{pnorm}) exists in R.
The argument \code{distribution} is the name of the cumulative distribution function with the \code{"p"} removed.

\code{zscoreGamma}, \code{zscoreT} and \code{zscoreHyper} are specific functions for the gamma, t and hypergeometric distributions respectively.

\code{tZscore} is the inverse of \code{zscoreT}, and computes t-distribution equivalents for standard normal deviates.

Care is taken to do the computations accurately in both tails of the distributions.
}

\author{Gordon Smyth}

\seealso{
\code{\link{qnorm}}, \code{\link{pgamma}}, \code{\link{pt}} in the stats package.
}

\examples{
# First three are equivalent
zscore(c(1,2.5), dist="gamma", shape=0.5, scale=2)
zscore(c(1,2.5), dist="chisq", df=1)
zscoreGamma(c(1,2.5), shape=0.5, scale=2)

zscoreT(2, df=3)
tZscore(2, df=3)
}

\keyword{distribution}
