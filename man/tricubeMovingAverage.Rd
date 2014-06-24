\name{tricubeMovingAverage}
\alias{tricubeMovingAverage}
\alias{tricubeMovingAverage}
\title{Moving Average Smoother With Tricube Weights}
\description{
Apply a moving average smoother with tricube distance weights to the columns of a matrix.
}
\usage{
tricubeMovingAverage(x,span=0.5,full.length=TRUE)
}
\arguments{
  \item{x}{numeric vector}
  \item{span}{proportion of points included in the local window}
  \item{full.length}{logical value, should output have same number of length as input?}
}
\details{
This function smooths a vector (considered as a time series) using a moving average with tricube weights.
This is similar to a loess curve of degree zero, with a couple of differences:
a continuity correction is applied when computing the neighbouring points and, when \code{full.length=TRUE}, the span halves at the end points.

The \code{filter} function in the stats package is called to do the low-level calculations.

This function is used by \code{\link{barcodeplot}} to compute enrichment worms.
}
\value{
Numeric vector of smoothed values.
If \code{full.length=TRUE}, of same length as \code{x}.
If \code{full.length=FALSE}, has \code{width-1} fewer rows than \code{x}.
}
\examples{
x <- rbinom(100,size=1,prob=0.5)
plot(1:100,tricubeMovingAverage(x))
}
\seealso{
\code{\link{filter}}, \code{\link{barcodeplot}}
}
\author{Gordon Smyth}
\keyword{smooth}
\concept{Loess curve}

