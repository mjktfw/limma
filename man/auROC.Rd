\name{auROC}
\alias{auROC}
\title{Area Under Receiver Operating Curve}
\description{
Compute exact area under the ROC for empirical data.
}
\usage{
auROC(truth, stat)
}
\arguments{
  \item{truth}{numeric vector of 0 and 1 indicating whether the null or alternative respectively is true for each case.
  If \code{stat} is missing then \code{truth} is assuming to be already sorted in decreasing test statistic order.}
  \item{stat}{numeric vector containing test statistics.
  Hypotheses are to be reject if \code{stat} exceeds a given threshold.}
}
\details{
The values in each column of \code{tstat} are assumed to follow a mixture of an ordinary t-distribution, with mixing proportion \code{1-proportion}, and \code{(v0+v1)/v1} times a t-distribution, with mixing proportion \code{proportion}.
}
\value{
Numeric vector giving area under the curve, 1 being perfect and 0 being the minimum.
\code{NULL} is truth has zero length.
}
\seealso{
\code{\link[ROC]{AUC}}
}
\examples{
auROC(c(1,1,0,0,0))
truth <- rbinom(30,size=1,prob=0.2)
stat <- rchisq(30,df=2)
auROC(truth,stat)
}
\author{Gordon Smyth}
\keyword{htest}
