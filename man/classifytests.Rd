\name{classifyTests}
\alias{classifyTests}
\alias{classifyTests.matrix}
\alias{classifyTests.vector}
\title{Treat Simultaneous T-Tests as Classification Problem}
\description{
Classify a series of related t-statistics as up, down or not significant.
}
\usage{
classifyTests(tstat,cor.matrix=diag(ntests),df.residual=Inf,p.value=0.01)
classifyTests.matrix(tstat,cor.matrix=diag(ntests),df.residual=Inf,p.value=0.01)
classifyTests.vector(tstat, Q=diag(p), pF=ncol(Q)*4)
}
\arguments{
  \item{tstat}{numeric matrix (\code{classifyTests}) or vector (\code{classifyTests.vector}) of t-statistics. Missing values are allowed.}
  \item{cor.matrix}{correlation matrix between the coefficients and the covariance matrix between the t-statistics}
  \item{df.residual}{numeric vector (\code{classifyTests}) or scalar (\code{classifyTests.vector}) giving the degrees of freedom for the t-statistics}
  \item{p.value}{numeric scalar between 0 and 1 giving the desired size of the test}
  \item{Q}{matrix spanning range space of \code{cor.matrix}}
  \item{pF}{critical value for F-type statistic}
}
\value{
  A matrix (\code{classifyTests}) or vector (\code{classifyTests.vector}) with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive.
}
\details{
F-statistics are constructed from the t-statistics.
If the overall F-statistic is significant, then the function makes a best choice as to which t-statistics contributed to this results.
The t-statistics are ranked.
Any t-statistic is called significant if the F-test is still significant for that row when all the larger t-statistics are set to the same absolute size as the t-statistic in question.

\code{classifyTests.vector} classifies a single vector of t-statistics, and is not intended to be called directly by the user.
\code{classifyTests} classifies each row of t-statistics.
}
\seealso{
\code{\link{ebayes}}
}
\author{Gordon Smyth}
\examples{
tstat <- matrix(rt(20*3,df=10),20,3)
tstat[1:10,1:2] <- tstat[1:10,1:2] + 4
classifyTests(tstat,diag(3),10,0.1)
}
\keyword{htest}
