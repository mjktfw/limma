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
\section{Methods}{
Methods are defined in limma for \code{object} an object of class
\code{\link[limma:normalizebetweenarrays-methods]{matrix}},
\code{\link[limma:normalizebetweenarrays-methods]{list}} or
\code{\link[limma:normalizebetweenarrays-methods]{MAList}}.

Type \code{showMethods("normalizeBetweenArrays")} at the R prompt for a complete list of classes for which methods are defined.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for normalization is given in \link{4.Normalization}.

See also \code{\link[affy:normalize-methods]{normalize}} in the affy package.
}
\keyword{multivariate}
