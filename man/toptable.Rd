\name{toptable}
\alias{toptable}
\title{Table of Top Genes from Linear Model Fit}
\description{
Extract a table of the top-ranked genes from a linear model fit.
}
\usage{
toptable(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="holm",sort.by="B",...)
}
\arguments{
  \item{fit}{output list from \code{lm.series}, \code{gls.series} or \code{rlm.series}}
  \item{coef}{column number of the effect or contrast to rank the genes on}
  \item{number}{how many genes to pick out}
  \item{genelist}{a data frame containing the gene allocation list or a vector containing the gene names}
  \item{A}{matrix of A-values or vector of average A-values.}
  \item{eb}{output list from \code{ebayes(fit)}}
  \item{adjust.method}{method to use to adjust the P-values for multiple testing, e.g., "holm" or "fdr". See \code{\link[ctest]{p.adjust}} for the available options. If \code{NULL} or \code{"none"} then the P-values are not adjusted.}
  \item{sort.by}{statistic to rank genes by.  Possibilities are \code{"M"}, \code{"A"}, \code{"T"}, \code{"P"} or \code{"B"}.}
  \item{...}{any other arguments are passed to \code{ebayes} if \code{eb} is \code{NULL}}
}
\value{
  A dataframe with a row for the \code{number} top genes and the following columns:
  \item{genelist}{if genelist was included as input}
  \item{M}{estimate of the effect or the contrast, on the log2 scale}
  \item{t}{moderated t-statistic}
  \item{P.Value}{nominal P-value}
  \item{B}{log odds that the gene is differentially expressed}
}
\details{
This function summarizes a fit object produced by \code{lm.series}, \code{gls.series} or \code{rlm.series} by selecting the top-ranked genes for any given contrast.
}
\seealso{
\code{\link{ebayes}}, \code{\link[base]{p.adjust}}, \code{\link{lm.series}}, \code{\link{gls.series}}, \code{\link{rlm.series}}.
}
\author{Gordon Smyth}
\examples{
#  Simulate gene expression data,
#  6 microarrays and 100 genes with first gene differentially expressed
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,1:3] <- M[1,1:3] + 2
#  Design matrix includes two treatments, one for first 3 and one for last 3 arrays
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- lm.series(M,design=design)
toptable(fit)
}
\keyword{htest}
