\name{lmscFit}
\alias{lmscFit}
\title{Single Channel Linear Model for Series of Arrays}
\description{Fit single channel linear model for each gene given a series of arrays}
\usage{
lmscFit(object, design, correlation)
}
\arguments{
  \item{object}{an \code{\link[limma:malist]{MAList}} object or a list from which \code{M} and \code{A} values may be extracted}
  \item{design}{a numeric matrix containing the design matrix for linear model in terms of the individual channels.
  The number of rows should be twice the number of arrays.
  The number of columns will determine the number of coefficients estimated for each gene.}
  \item{correlation}{numeric value giving the intra-spot correlation}
}
\value{
An object of class \code{\link[limma:marraylm]{MArrayLM}}
}
\details{
For two color arrays, the channels measured on the same set of arrays are correlated.
The \code{M} and \code{A} however are uncorrelated for each gene.
This function fits a linear model to the set of M and A-values for each gene after re-scaling the M and A-values to have equal variances.
The input correlation determines the scaling required.

Missing values in \code{M} or \code{A} are not allowed.
}
\author{Gordon Smyth}
\seealso{
\code{\link[base:lmfit]{lm.fit}}.

An overview of methods for single channel analysis in limma is given by \link{6.SingleChannel}.
}
\keyword{models}
\keyword{regression}
