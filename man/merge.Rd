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
\section{Methods}{
A method is defined in limma for \code{x} and \code{y} objects of class \code{\link[limma:merge.rglist]{RGList}}.
The default method \code{\link[base]{merge.default}} defined in R base is designed for data.frames.

Type \code{showMethods("merge")} at the R prompt for a complete list of classes for which methods are defined.
}
\author{Gordon Smyth}
\keyword{models}
