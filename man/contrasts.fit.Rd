\name{contrasts.fit}
\alias{contrasts.fit}
\title{Contrast Information from Linear Model Fit}
\description{
Given an \code{lm.series} fit for a oneway model, compute estimated coefficients and standard errors for a given set of contrasts.
}
\usage{
contrasts.fit(fit,contrasts)
}
\arguments{
  \item{fit}{object produced by the function \code{lm.series} or equivalent. A list containing components \code{coefficients} and \code{stdev.unscaled}.}
  \item{contrasts}{matrix with columns containing contrasts. May be a vector if there is only one contrast.}
}
\value{
  An object of the same type as produced by \code{lm.series}. This is a list components
  \item{coefficients}{numeric matrix containing the estimated coefficients for each contrasts for each gene.}
  \item{stdev.unscaled}{numeric matrix conformal with \code{coef} containing the unscaled standard deviations for the coefficient estimators.}
  \item{...}{any other components input in \code{fit}}
}
\details{
This function accepts input from any of the functions \code{lm.series}, \code{rlm.series} or \code{gls.series}.
The design matrix used for this fit must have orthogonal columns.

The idea is to fit a saturated oneway model using of the above functions, then use \code{contrasts.fit} to obtain coefficients and standard errors for any number of contrasts of the coefficients of the oneway model.
}
\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\examples{
#  Simulate gene expression data,
#  6 microarrays and 100 genes with one gene differentially expressed in first 3 arrays
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,1:3] <- M[1,1:3] + 2
#  Design matrix corresponds to oneway layout, columns are orthogonal
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- lm.series(M,design=design)
#  Would like to consider original two estimates plus difference between first 3 and last 3 arrays
contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),"Last3-First3"=c(-1,1))
fit2 <- contrasts.fit(fit,contrasts=contrast.matrix)
eb <- ebayes(fit2)
#  Large values of eb$t indicate differential expression
clas <- classifyTests(eb$t,design=design,contrasts=contrast.matrix,df=fit2$df+eb$df)
}
\keyword{htest}
