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
\section{Methods}{
Methods are defined in limma for \code{object} an object of class
\code{\link[limma:normalizeWithinArrays-methods]{list}},
\code{\link[limma:normalizeWithinArrays-methods]{RGList}} or
\code{\link[limma:normalizeWithinArrays-methods]{MAList}}.

Type \code{showMethods("normalizeBetweenArrays")} at the R prompt for a complete list of classes for which methods are defined.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for normalization is given in \link{4.Normalization}.

See also \code{\link{normalizeBetweenArrays}} and \code{\link[marrayNorm]{maNorm}} in the marrayNorm package.
}
\keyword{models}
