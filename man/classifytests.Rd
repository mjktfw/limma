\name{classifyTests}
\alias{classifyTests}
\title{Treat Simultaneous T-Tests as Classification Problem}
\description{
Classify a series of related t-statistics as up, down or not significant.
}
\usage{
classifyTests(tstat, cor.matrix=NULL, design=NULL, contrasts=NULL, df=Inf, p.value=0.01)
}
\arguments{
  \item{tstat}{numeric matrix of t-statistics or an \code{MArrayLM} object from which the t-statistics may be extracted.}
  \item{cor.matrix}{covariance matrix of each row of t-statistics.  Defaults to the identity matrix.}
  \item{design}{full rank numeric design matrix.  Not used if \code{cor.matrix} is specified.}
  \item{contrasts}{numeric matrix with columns specifying contrasts of the coefficients of interest.  Not used if \code{cor.matrix} is specified.}
  \item{df}{numeric vector giving the degrees of freedom for the t-statistics.
  May have length 1 or length equal to the number of rows of \code{tstat}.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test}
}
\value{
  A matrix with elements \code{-1}, \code{0} or \code{1} depending on whether each t-statistic is classified as significantly negative, not significant or significantly positive respectively.
}
\details{
Rows of \code{tstat} correspond to genes and columns to coefficients or contrasts.
For each row of \code{tstat}, F-statistics are constructed from the t-statistics.
If the overall F-statistic is significant, then the function makes a best choice as to which t-statistics contributed to this result.
The methodology is based on the principle that any t-statistic should be called significant if the F-test is still significant for that row when all the larger t-statistics are set to the same absolute size as the t-statistic in question.

If \code{tstat} is an \code{MArrayLM} object, then all arguments except for \code{p.value} are extracted from it.

\code{cor.matrix} is the same as the correlation matrix of the coefficients from which the t-statistics are calculated.
If \code{cor.matrix} is not specified, then it is calculated from \code{design} and \code{contrasts} if at least \code{design} is specified or else defaults to the identity matrix.
In terms of \code{design} and \code{contrasts}, \code{cor.matrix} is obtained by standardizing the matrix
\code{ t(contrasts) \%*\% solve(t(design) \%*\% design) \%*\% contrasts }
to a correlation matrix.
}
\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\examples{
tstat <- matrix(c(0,5,0, 0,2.5,0, -2,-2,2, 1,1,1), 4, 3, byrow=TRUE)
classifyTests(tstat)

# See also the example for contrasts.fit
}
\keyword{htest}
