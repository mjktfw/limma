\name{tmixture}
\alias{tmixture.vector}
\alias{tmixture.matrix}
\title{Estimate Scale Factor in Mixture of t-Distributions}
\description{
This function estimates the unscaled standard deviation of the log fold change for differentially expressed genes.
It is called by the function \code{ebayes} and is not intended to be called by users.
}
\usage{
tmixture.vector(tstat,stdev.unscaled,df,proportion,c0lim=NULL)
tmixture.matrix(tstat,stdev.unscaled,df,proportion,c0lim=NULL)
}
\arguments{
  \item{tstat}{numeric vector or matrix of t-statistics}
  \item{stdev.unscaled}{numeric matrix conformal with \code{tstatf} containing the unscaled standard deviations for the coefficient estimators}
  \item{df}{numeric vector giving the degrees of freedom associated with \code{tstat}}
  \item{proportion}{assumed proportion of genes which are differentially expressed}
  \item{c0lim}{optional upper limit for the estimated unscaled standard deviation}
}
\value{
Numeric vector of length equal to the number of columns of \code{tstat} and \code{stdev.unscaled}.
}
\seealso{
\code{\link{ebayes}}
}
\author{Gordon Smyth}
\keyword{htest}
