\name{ebayes}
\alias{ebayes}
\alias{eBayes}
\alias{treat}
\title{Empirical Bayes Statistics for Differential Expression}
\description{Given a microarray linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a common value.}
\usage{
ebayes(fit, proportion=0.01, stdev.coef.lim=c(0.1,4), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
eBayes(fit, proportion=0.01, stdev.coef.lim=c(0.1,4), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
treat(fit, lfc=0, trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit}.
  For \code{ebayes} only, \code{fit} can alternatively be an unclassed list produced by \code{lm.series}, \code{gls.series} or \code{mrlm} containing components \code{coefficients}, \code{stdev.unscaled}, \code{sigma} and \code{df.residual}.}
  \item{proportion}{numeric value between 0 and 1, assumed proportion of genes which are differentially expressed}
  \item{stdev.coef.lim}{numeric vector of length 2, assumed lower and upper limits for the standard deviation of log2-fold-changes for differentially expressed genes}
  \item{trend}{logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.}
  \item{robust}{logical, should the estimation of \code{df.prior} and \code{var.prior} be robustified against outlier sample variances?}
  \item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail proportions of \code{x} to Winsorize. Used only when \code{robust=TRUE}.}
  \item{lfc}{the minimum log2-fold-change that is considered scientifically meaningful}
}
\value{
\code{eBayes} produces an object of class \code{MArrayLM} (see \code{\link{MArrayLM-class}}) containing everything found in \code{fit} plus the following added components:
  \item{t}{numeric vector or matrix of moderated t-statistics}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{s2.prior}{estimated prior value for \code{sigma^2}. A vector if \code{covariate} is non-\code{NULL}, otherwise a scalar.}
  \item{df.prior}{degrees of freedom associated with \code{s2.prior}}
  \item{df.total}{numeric vector of total degrees of freedom associated with t-statistics and p-values. Equal to \code{df.prior+df.residual} or \code{sum(df.residual)}, whichever is smaller.}
  \item{s2.post}{numeric vector giving the posterior values for \code{sigma^2}}
  \item{lods}{numeric vector or matrix giving the log-odds of differential expression}
  \item{var.prior}{estimated prior value for the variance of the log2-fold-change for differentially expressed gene}
  \item{F}{numeric vector of moderated F-statistics for testing all contrasts defined by the columns of \code{fit} simultaneously equal to zero}
  \item{F.p.value}{numeric vector giving p-values corresponding to \code{F}}
\code{treat} a produces an \code{MArrayLM} object similar to \code{eBayes} but without \code{lods}, \code{var.prior}, \code{F} or \code{F.p.value}.
\code{ebayes} produces an ordinary list containing the above components except for \code{F} and \code{F.p.value}.
}
\details{
These functions is used to rank genes in order of evidence for differential expression.
They use an empirical Bayes method to shrink the probe-wise sample variances towards a common value and to augmenting the degrees of freedom for the individual variances (Smyth, 2004).
The functions accept as input argument \code{fit} a fitted model object from the functions \code{lmFit}, \code{lm.series}, \code{mrlm} or \code{gls.series}.
The fitted model object may have been processed by \code{contrasts.fit} before being passed to \code{eBayes} to convert the coefficients of the design matrix into an arbitrary number of contrasts which are to be tested equal to zero.
The columns of \code{fit} define a set of contrasts which are to be tested equal to zero.

The empirical Bayes moderated t-statistics test each individual contrast equal to zero.
For each probe (row), the moderated F-statistic tests whether all the contrasts are zero.
The F-statistic is an overall test computed from the set of t-statistics for that probe.
This is exactly analogous the relationship between t-tests and F-statistics in conventional anova, except that the residual mean squares and residual degrees of freedom have been moderated between probes.

The estimates \code{s2.prior} and \code{df.prior} are computed by \code{fitFDist}.
\code{s2.post} is the weighted average of \code{s2.prior} and \code{sigma^2} with weights proportional to \code{df.prior} and \code{df.residual} respectively.
The \code{lods} is sometimes known as the B-statistic.
The F-statistics \code{F} are computed by \code{classifyTestsF} with \code{fstat.only=TRUE}.

\code{eBayes} doesn't compute ordinary (unmoderated) t-statistics by default, but these can be easily extracted from 
the linear model output, see the example below.

\code{ebayes} is the earlier and leaner function, kept for background capatability, while
\code{eBayes} is the later more object-orientated version.
The difference is that \code{ebayes} outputs only the empirical Bayes statistics whereas \code{eBayes} adds them to the fitted model object \code{fit}.
\code{eBayes} is recommended for routine use as it produces objects containing all the necessary components for downstream analysis

\code{treat} computes empirical Bayes moderated-t p-values relative to a minimum required fold-change threshold.
Use \code{\link{topTreat}} to summarize output from \code{treat}.
Instead of testing for genes which have log-fold-changes different from zero, it tests whether the log2-fold-change is greater than \code{lfc} in absolute value (McCarthy and Smyth, 2009).
\code{treat} is concerned with p-values rather than posterior odds, so it does not compute the B-statistic \code{lods}.
The idea of thresholding doesn't apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed.

If \code{trend=TRUE} then an intensity-dependent trend is fitted to the prior variances \code{s2.prior}.
Specifically, \code{squeezeVar} is called with the \code{covariate} equal to \code{Amean}, the average log2-intensity for each gene.
See \code{\link{squeezeVar}} for more details.

If \code{robust=TRUE} then the robust empirical Bayes procedure of Phipson et al (2013) is used.
See \code{\link{squeezeVar}} for more details.
}
\seealso{
\code{\link{squeezeVar}}, \code{\link{fitFDist}}, \code{\link{tmixture.matrix}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth and Davis McCarthy}
\references{
McCarthy, D. J., and Smyth, G. K. (2009). Testing significance relative to a fold-change threshold is a TREAT. \emph{Bioinformatics}.
\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp053}

Loennstedt, I., and Speed, T. P. (2002). Replicated microarray data. \emph{Statistica Sinica} \bold{12}, 31-46.

Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2013).
Empirical Bayes in the presence of exceptional cases, with application to microarray data.
Bioinformatics Division, Walter and Eliza Hall Institute of Medical Research, Melbourne, Australia.
\url{http://www.statsci.org/smyth/pubs/RobustEBayesPreprint.pdf}

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume \bold{3}, Article 3.
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
}
\examples{
#  See also lmFit examples

#  Simulate gene expression data,
#  6 microarrays and 100 genes with one gene differentially expressed
set.seed(2004); invisible(runif(100))
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,] <- M[1,] + 1
fit <- lmFit(M)

#  Moderated t-statistic
fit <- eBayes(fit)
topTable(fit)

#  Ordinary t-statistic
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma

#  Q-Q plots of t statistics
#  Points off the line may be differentially expressed
par(mfrow=c(1,2))
qqt(ordinary.t, df=fit$df.residual, main="Ordinary t")
abline(0,1)
qqt(fit$t, df=fit$df.total,main="Moderated t")
abline(0,1)
par(mfrow=c(1,1))
}
\keyword{htest}
