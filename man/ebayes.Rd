\name{ebayes}
\alias{ebayes}
\alias{eBayes}
\title{Empirical Bayes Statistics for Differential Expression}
\description{Given a series of related parameter estimates and standard errors, compute moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.}
\usage{
ebayes(fit,proportion=0.01,stdev.coef.lim=c(0.1,4))
eBayes(fit,proportion=0.01,stdev.coef.lim=c(0.1,4))
}
\arguments{
  \item{fit}{a list object produced by \code{lm.series}, \code{gls.series}, \code{rlm.series} or \code{lmFit} containing components \code{coefficients}, \code{stdev.unscaled}, \code{sigma} and \code{df.residual}}
  \item{proportion}{numeric value between 0 and 1, assumed proportion of genes which are differentially expressed}
  \item{stdev.coef.lim}{numeric vector of length 2, assumed lower and upper limits for the standard deviation of log2 fold changes for differentially expressed genes}
}
\value{
\code{ebayes} produces an ordinary list with the following components.
\code{eBayes} adds the following components to \code{fit} to produce an augmented object, usually of class \code{MArrayLM}.
  \item{t}{numeric vector or matrix of penalized t-statistics}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{s2.prior}{estimated prior value for \code{sigma^2}}
  \item{df.prior}{degrees of freedom associated with \code{s2.prior}}
  \item{s2.post}{vector giving the posterior values for \code{sigma^2}}
  \item{lods}{numeric vector or matrix giving the log-odds of differential expression}
  \item{var.prior}{estimated prior value for the variance of the log2-fold-change for differentially expressed gene}
}
\details{This function is used to rank genes in order of evidence for differential expression.
The function accepts as input output from the functions \code{lm.series}, \code{rlm.series} or \code{gls.series}.
The estimates \code{s2.prior} and \code{df.prior} are computed by \code{fdist.fit}.
\code{s2.post} is the weighted average of \code{s2.prior} and \code{sigma^2} with weights proportional to \code{df.prior} and \code{df.residual} respectively.

The \code{lods} is sometimes known as the B-statistic.
}
\seealso{
\code{\link{fitFDist}}, \code{\link{tmixture.matrix}}.

An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth}
\references{
Lönnstedt, I. and Speed, T. P. (2002). Replicated microarray data. \emph{Statistica Sinica} \bold{12}, 31-46.

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3}, No. 1, Article 3. \url{http://www.bepress.com/sagmb/vol3/iss1/art3}
}
\examples{
#  Simulate gene expression data,
#  6 microarrays and 100 genes with one gene differentially expressed
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,] <- M[1,] + 1.6
fit <- lm.series(M)
eb <- ebayes(fit)
qqt(eb$t,df=eb$df+fit$df)
abline(0,1)
#  Points off the line may be differentially expressed
}
\keyword{htest}
