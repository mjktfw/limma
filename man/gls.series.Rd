\name{gls.series}
\alias{gls.series}
\title{Generalized Least Squares for Series of Microarrays}
\description{Fit linear models for each gene to a series of microarrays. Fit is by generalized least squares allowing for correlation between duplicate spots.}
\usage{gls.series(M,design=rep(1,ncol(M)),ndups=2,spacing=1,correlation=NULL,weights=NULL,...)}
\arguments{
  \item{M}{a numeric matrix. Usually the log-ratios of expression for a series of cDNA microarrrays with rows corresponding to genes and columns to arrays.}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{M}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array. \code{nrow(M)} must be divisible by \code{ndups}.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{correlation}{the inter-duplicate correlation.}
  \item{weights}{an optional numeric matrix of the same dimension as \code{M} containing weights for each spot. If it is of different dimension to \code{M}, it will be filled out to the same size.}
  \item{...}{other optional arguments to be passed to \code{dupcor.series}.}
}
\value{
  A list with components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each linear model. Same number of rows as \code{M}, same number of columns as \code{design}.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators. The standard errors are given by \code{stdev.unscaled * sigma}.}
  \item{sigma}{numeric vector containing the residual standard deviation for each gene.}
  \item{df.residual}{numeric vector giving the degrees of freedom corresponding to \code{sigma}.}
  \item{correlation}{inter-duplicate correlation.}
}
\details{
Normally \code{dupcor.series} will be called before \code{gls.series} to estimate the inter-duplicate correlation. Then t-statistics will be formed to rank the genes in order of evidence for differential expression.
}
\seealso{
\code{\link{dupcor.series}}.

An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\examples{
M <- matrix(rnorm(10*6),10,6)
cor.out <- dupcor.series(M)
#  See cor.out$cor for estimated common correlation
fit <- gls.series(M,correlation=cor.out$cor)
}
\keyword{models}
\keyword{regression}

