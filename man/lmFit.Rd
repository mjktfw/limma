\name{lmFit}
\alias{lmFit}
\title{Linear Model for Series of Arrays}
\description{Fit linear model for each gene given a series of arrays}
\usage{
lmFit(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation=0.75,weights=NULL,method="ls",...) 
}
\arguments{
  \item{object}{object of class \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm}, \code{exprSet} or \code{PLMset} containing log-ratios or log-values of expression for a series of microarrays}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to coefficients to be estimated.  Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{positive integer giving the number of times each gene is printed on an array}
  \item{spacing}{positive integer giving the spacing between duplicate spots, \code{spacing=1} for consecutive spots}
  \item{block}{vector or factor specifying a blocking variable on the arrays. Has length equal to the number of arrays.}
  \item{correlation}{the inter-duplicate or inter-technical replicate correlation}
  \item{weights}{optional numeric matrix containing weights for each spot}
  \item{method}{character string, \code{"ls"} for least squares or \code{"robust"} for robust regression}
  \item{...}{other optional arguments to be passed to \code{lm.series}, \code{gls.series} or \code{rlm.series}}
}

\value{
Object of class \code{\link[limma:marraylm]{MArrayLM}}
}

\details{
A linear model is fitted for each gene by calling one of the lower-level functions \code{lm.series}, \code{gls.series} or \code{rlm.series}.
The arguments \code{design}, \code{ndups}, \code{spacing} and \code{weights} will be extracted from the data \code{object} if available and do not normally need to set explicitly in the call; if any of these are set in the call then they will over-ride the slots or components in the data \code{object}.

If \code{object} is an \code{PLMset}, then expression values will be taken from the slot \code{chip.coefs} and weights will be computed from \code{se.chip.coefs}.
If \code{object} is an \code{exprSet}, then expression values will be taken from the \code{exprs} slot, but weights will not be computed.
}

\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}

\author{Gordon Smyth}
\keyword{models}
\keyword{regression}
