\name{dupcor.series}
\alias{duplicateCorrelation}
\alias{dupcor.series}
\title{Correlation Between Duplicates}
\description{Estimate the correlation between duplicate spots (replicate spots on the same array) from a series of arrays.}
\usage{
duplicateCorrelation(object,design=rep(1,ncol(M)),ndups=2,spacing=1,trim=0.15,weights=NULL)
dupcor.series(M,design=rep(1,ncol(M)),ndups=2,spacing=1,initial=0.8,trim=0.15,weights=NULL)
}
\arguments{
  \item{object}{a numeric matrix of log-ratios or an \code{\link[limma:malist]{MAList}} object from which the log-ratios can be extracted.
  If \code{object} is an \code{MAList} then the arguments \code{design}, \code{ndups}, \code{spacing} and \code{weights} will be extracted from it if available and do not have to be specified as arguments.}
  \item{M}{a numeric matrix. Usually the log-ratios of expression for a series of cDNA microarrrays with rows corresponding to genes and columns to arrays.}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{M}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array. \code{nrow(M)} must be divisible by \code{ndups}.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{initial}{a numeric value between -1 and 1 giving an initial estimate for the correlation.
  Not currently used.}
  \item{trim}{the fraction of observations to be trimmed from each end of \code{tanh(cor.genes)} when computing the trimmed mean.}
  \item{weights}{an optional numeric matrix of the same dimension as \code{M} containing weights for each spot. If smaller than \code{M} then it will be filled out the same size.}
}
\value{
  A list with components
  \item{cor}{the average estimated inter-duplicate correlation. The average is the 0.1 trimmed mean of the correlations for individual genes on the tanh-transformed scale.}
  \item{cor.genes}{a numeric vector of length \code{nrow(M)/ndups} giving the individual gene correlations.}
}
\details{
This function estimates the between-duplicate correlation using REML individually for each gene.
It also returns a robust average of the individual correlations which can be used as input for 
functions \code{lmFit} or \code{gls.series}.

\code{dupcor.series} is a depreciated function which may be removed in a future version.
It produces the same value as \code{duplicateCorrelation}.
}
\note{
This function may take long time to execute as it makes a call to \code{\link[nlme]{gls}} for each gene.
Execution could be speeded up greatly if it could be assumed that \code{M} contains no NAs.
}
\seealso{
These functions use \code{\link[statmod]{randomizedBlockFit}} in the statmod package.

An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\references{
Smyth, G. K., Michaud, J., and Scott, H. (2003). The use of within-array duplicate spots for assessing differential expression in microarray experiments.
\url{http://www.statsci.org/smyth/pubs/dupcor.pdf}
}
\examples{
#  See gls.series for an example
}
\keyword{multivariate}
