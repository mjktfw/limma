\name{merge-generic}
\alias{merge}
\title{Merge Data Objects - generic}
\description{
A generic function which merges data objects to produce a single object with the combined data.
}
\usage{
merge(x,y,...)
}
\arguments{
  \item{x}{a data object}
  \item{y}{a data object of the same class as \code{x} to be merged with \code{x}}
  \item{...}{any other arguments}
}
\author{Gordon Smyth}
\seealso{
The default method \code{\link[base]{merge.default}} defined in R base is designed for data.frames.

\code{\link{merge.RGList.RGList}}.
Type \code{showMethods("merge")} at the R prompt to see what other methods are available.
}
\keyword{models}
