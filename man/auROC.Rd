\name{auROC}
\alias{auROC}
\title{Area Under Receiver Operating Curve}
\description{
Compute exact area under the ROC for empirical data.
}
\usage{
auROC(truth, stat=NULL)
}
\arguments{
  \item{truth}{logical vector, or numeric vector of 0s and 1s, indicating whether each case is a true positive.}
  \item{stat}{numeric vector containing test statistics used to rank cases, from largest to smallest.
  If \code{NULL}, then \code{truth} is assumed to be already sorted in decreasing test statistic order.}
}
\details{
A receiver operating curve (ROC) is a plot of sensitivity (true positive rate) versus error (false positive rate) for a statistical test or binary classifier.
The area under the ROC is a well accepted measure of test performance.
It is equivalent to the probability that a randomly chosen pair of cases is corrected ranked.

Here we consider a test statistic \code{stat}, with larger values being more significant, and a vector \code{truth} indicating whether the null hypothesis is in fact true.
Correct ranking here means that \code{truth[i]} is greater than or equal to \code{truth[j]} when \code{stat[i]} is greater than \code{stat[j]}. 
The function computes the exact area under the empirical ROC curve defined by \code{truth} when ordered by \code{stat}.
}
\value{
Numeric vector giving area under the curve, 1 being perfect and 0 being the minimum, or \code{NULL} if \code{truth} has zero length.
}
\seealso{
See \link{08.Tests} for other functions for testing and processing p-values.
}
\examples{
auROC(c(1,1,0,0,0))
truth <- rbinom(30,size=1,prob=0.2)
stat <- rchisq(30,df=2)
auROC(truth,stat)
}
\author{Gordon Smyth}
\keyword{htest}
