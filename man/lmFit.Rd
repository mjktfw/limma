\name{lmFit}
\alias{lmFit}
\title{Linear Model for Series of Arrays}
\description{Fit linear model for each gene given a series of arrays}
\usage{
lmFit(object,design=NULL,ndups=1,spacing=1,correlation=0.75,weights=NULL,method="ls",...) 
}
\arguments{
  \item{object}{object of class \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm} or \code{exprSet} containing log-ratios or log-values of expression for a series of microarrays}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to coefficients to be estimated.  Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array}
  \item{spacing}{the spacing between duplicate spots, \code{spacing=1} for consecutive spots}
  \item{correlation}{the inter-duplicate correlation}
  \item{weights}{an optional numeric matrix containing weights for each spot}
  \item{method}{character string, \code{"ls"} for least squares or \code{"robust"} for robust regression}
  \item{...}{other optional arguments to be passed to \code{lm.series}, \code{gls.series} or \code{rlm.series}}
}

\value{
Object of class \code{\link[limma:marraylm]{MArrayLM}}
}

\details{
A linear model is fitted for each gene by calling one of \code{lm.series}, \code{gls.series} or \code{rlm.series}.
Note that the arguments \code{design}, \code{ndups} will be extracted from the data \code{object} if available and do not normally need to set explicitly in the call.
If arguments are set in the call then they will over-ride slots or components in the data \code{object}.
}

\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}

\author{Gordon Smyth}
\keyword{models}
\keyword{regression}
