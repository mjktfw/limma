\name{normalizeWithinArrays}
\alias{normalizeWithinArrays}
\title{Normalize Within Arrays - generic}
\description{
A generic function which normalizes microarray data.
Normalization is intended to remove from the expression measures any systematic trends which arise from the microarray technology rather than from differences between the probes or between the target RNA samples hybridized to the arrays.
}
\usage{
normalizeWithinArrays(object,...)
}
\arguments{
  \item{object}{a data object containing microarray data}
  \item{...}{any other arguments}
}
\author{Gordon Smyth}
\seealso{
Type \code{showMethods("normalizeWithinArrays")} at the R prompt to see what methods are available.
Help on individual methods is generally available as \code{normalizeWithinArrays}.<class> where <class> is the class of the data object.

See also \code{\link{normalizeBetweenArrays}}, \code{\link[marrayNorm]{maNorm}} in the marrayNorm package and \code{\link[affy]{normalize}} in the affy package.
}
\keyword{models}
