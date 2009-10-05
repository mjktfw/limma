\name{zscore}
\alias{zscoreGamma}
\alias{zscoreT}
\alias{tZscore}

\title{Z-score Equivalents}

\description{
Compute z-score equivalents of gamma or t-distribution random deviates.
}

\usage{
zscoreGamma(q, shape, rate = 1, scale = 1/rate) 
zscoreT(x, df)
tZscore(x, df)
}

\arguments{
\item{q, x}{numeric matrix for vector giving deviates of a random variaable}
\item{shape}{gamma shape parameter (>0)}
\item{rate}{gamma rate parameter (>0)}
\item{scale}{gamma scale parameter (>0)}
\item{df}{degrees of freedom (>0 for \code{zscoreT} or >=1 for \code{tZscore})}
}

\value{
Numeric vector giving equivalent deviates from a standard normal distribution (\code{zscoreGamma} and \code{zscoreT}) or deviates from a t-distribution (\code{tZscore}).
}

\details{
These functions compute the standard normal deviates which have the same quantiles as the given values in the specified distribution.
For example, if \code{z <- zscoreT(x,df=df)} then \code{pnorm(z)} equals \code{pt(x,df=df)}.
\code{tZscore} is the inverse of \code{zscoreT}.

Care is taken to do the computations accurately in both tails of the distributions.
}

\author{Gordon Smyth}
\seealso{
\code{\link{qnorm}}, \code{\link{pgamma}}, \code{\link{pt}} in the standard distribution packages.
}
\examples{
zscoreGamma(1, shape=1, scale=1)
zscoreT(2, df=3)
tZscore(2, df=3)
}
\keyword{distribution}
