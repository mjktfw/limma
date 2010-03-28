\name{avereps}
\alias{avereps}
\alias{avereps.default}
\alias{avereps.MAList}
\title{Average Over Irregular Replicate Spots}
\description{
Condense a microarray data object so that values for within-array replicate spots are replaced with their average.
}
\usage{
\method{avereps}{default}(x, ID=rownames(x))
\method{avereps}{MAList}(x, ID=NULL)
}
\arguments{
  \item{x}{a matrix-like object, usually a matrix or an \code{MAList} object.}
  \item{ID}{probe identifier.}
}
\details{
A new data object is computed in which each probe is represented by the average of its replicate spots.
For an \code{MAList} object, the components \code{M} and \code{A} are both averaged in this way, as \code{weights} and any matrices found in \code{object$other}.

For an \code{MAList} object, \code{ID} defaults to \code{MA$genes$ID} is that exists, otherwise to \code{rownames(MA$M)}.

If \code{x} is of mode \code{"character"}, then the replicate values are assumed to be equal and the first is taken as the average.
}
\value{
A data object of the same class as \code{x} with a row for each unique value of \code{ID}.
}
\author{Gordon Smyth}
\seealso{
  \code{\link{avedups}}. Also \code{\link[base]{rowsum}} in the base package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}

\examples{
x <- matrix(rnorm(8*3),8,3)
colnames(x) <- c("S1","S2","S3")
rownames(x) <- c("b","a","a","c","c","b","b","b")
avereps(x)
}
\keyword{array}
