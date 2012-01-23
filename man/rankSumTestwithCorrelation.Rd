\name{rankSumTestwithCorrelation}
\alias{rankSumTestwithCorrelation}

\title{Two Sample Wilcoxon-Mann-Whitney Rank Sum Test Allowing For Correlation}

\description{
A extension of the well-known rank-based test, but allowing for correlations between cases.
}

\usage{
rankSumTestwithCorrelation(index, statistics, correlation=NULL, vif=NULL, df=Inf)
}

\arguments{
  \item{index}{any index vector such that \code{statistics[index]} contains the values of the statistic for the test group.}
  \item{statistics}{numeric vector giving values of the test statistic.}
  \item{correlation}{numeric scalar, average correlation between cases in the test group.  Remaining cases are assumed independent.} 
  \item{vif}{variance inflation factor for the mean, equal to \code{1+(m-1)*correlation} where \code{m} is the number of genes in the test group.} 
  \item{df}{degrees of freedom which the correlation has been estimated.} 
}

\value{
Numeric vector containing \code{left.tail} and \code{right.tail} p-values respectively.
}

\details{
Tests whether the mean rank of statistics in the test group is greater or less than the mean rank of the remaining statistic values.

When the correlation (or variance inflation factor) is zero, the function performs the usual two-sample Wilcoxon-Mann-Whitney test rank sum test.
The test is implemented following the formulas given in Zar (1999) Section 8.10, including corrections for ties and for continuity.

The test allows for the possibility that cases in the test group may be more highly correlated on average than cases not in the group.
When the correlation is non-zero, the variance of the z-statistic is inflated, and the standard normal test is replaced with a t-distribution test.
When the correlation is positive, the VIF is positive and test will be more conservative.

This correlation-adjusted test was proposed and described by Wu and Smyth (2012).
}

\seealso{
\code{\link{wilcox.test}} performs the usual Wilcoxon-Mann-Whitney test assuming independence.

An overview of tests in limma is given in \link{08.Tests}.
}

\author{Gordon Smyth and Di Wu}

\references{
Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
Submitted.

Zar, JH (1999). \emph{Biostatistical Analysis 4th Edition}.
Prentice-Hall International, Upper Saddle River, New Jersey.
}

\examples{
stat <- rnorm(100)
index <- 1:10
stat[index] <- stat[1:10]+1

rankSumTestwithCorrelation(index, stat)
rankSumTestwithCorrelation(index, stat, correlation=0.1)

group <- rep(1,100)
group[index] <- 2
group <- factor(group)
wilcox.test(stat ~ group)
}

\keyword{htest}

