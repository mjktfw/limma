\name{normalizeBetweenArrays-generic}
\alias{normalizeBetweenArrays}
\title{Normalize Between Arrays - generic}
\description{
A generic function which normalizes expression intensities across a series of arrays.
}
\usage{
normalizeBetweenArrays(object,...)
}
\arguments{
  \item{object}{an object containing expression ratios for a series of arrays}
  \item{...}{any other arguments}
}
\author{Gordon Smyth}
\seealso{
\code{\link{normalizeBetweenArrays.matrix}}, \code{\link{normalizeBetweenArrays.MAList}}.
Type \code{showMethods("normalizeBetweenArrays")} at the R prompt to see what other methods are available.
}
\keyword{multivariate}
