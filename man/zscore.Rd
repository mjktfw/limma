\name{zscore}
\alias{zscoreGamma}
\alias{zscoreT}

\title{z-score equivalents}

\description{
Compute z-score equivalents of for gamma or t-distribution random deviates.
}

\usage{
zscoreGamma(q, shape, rate = 1, scale = 1/rate) 
zscoreT(x, df)
}

\arguments{
\item{q, x}{numeric matrix for vector giving deviates of a random variaable}
\item{shape}{gamma shape parameter}
\item{rate}{gamma rate parameter}
\item{scale}{gamma scale parameter}
\item{df}{degrees of freedom}
}

\value{
Numeric vector giving equivalent deviates from a standard normal distribution.
}

\details{
This function computes the deviates of a standard normalization distribution which have the same quantiles as the given values in the specified distribution.
For example, if \code{z <- zscoreT(x,df=df)} then \code{pnorm(z)} equals \code{pt(x,df=df)}.

Care is taken to do the computations accurately in both tails of the distributions.
}

\author{Gordon Smyth}
\seealso{
\code{\link[base]{qnorm}}, \code{\link[base]{pgamma}}, \code{\link[base]{pt}}
}
\examples{
zscoreGamma(1, shape=1, scale=1)
zscoreT(0, df=3)
}
\keyword{distribution}
