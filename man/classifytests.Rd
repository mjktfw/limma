\name{classifyTests}
\alias{classifyTests}
\title{Treat Simultaneous T-Tests as Classification Problem}
\description{
Classify a series of related t-statistics as up, down or not significant.
}
\usage{
classifyTests(tstat, cor.matrix=diag(ntests), df.residual=Inf, p.value=0.01)
}
\arguments{
  \item{tstat}{numeric matrix of t-statistics. Missing values are allowed.}
  \item{cor.matrix}{covariance matrix of the t-statistics.  Same as the correlation matrix of the corresponding coefficients.}
  \item{df.residual}{numeric vector (\code{classifyTests}) giving the degrees of freedom for the t-statistics}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
}
\value{
  A matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.
}
\details{
Rows of \code{tstat} correspond to genes and columns to coefficients or contrasts.
For each row of \code{tstat}, F-statistics are constructed from the t-statistics.
If the overall F-statistic is significant, then the function makes a best choice as to which t-statistics contributed to this result.
The t-statistics are ranked.
Any t-statistic is called significant if the F-test is still significant for that row when all the larger t-statistics are set to the same absolute size as the t-statistic in question.
}
\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\examples{
tstat <- matrix(c(0,5,0, 0,2.5,0, -2,-2,2, 1,1,1), 4, 3, byrow=TRUE)
classifyTests(tstat)
}
\keyword{htest}
