\name{squeezeVar}
\alias{squeezeVar}

\title{Squeeze Sample Variances}
\description{
Squeeze a set of sample variances together by computing empirical Bayes posterior means.
}

\usage{
squeezeVar(var, df, covariate=NULL)
}

\arguments{
\item{var}{numeric vector of independent sample variances.}
\item{df}{numeric vector of degrees of freedom for the sample variances.}
\item{covariate}{if non-\code{NULL}, \code{var.prior} will depend on this numeric covariate. Otherwise, \code{var.prior} is constant.}
}

\details{
The sample variances \code{var} are assumed to follow scaled chi-squared distributions.
An inverse chi-squared prior is assumed for the true variances.
The scale and degrees of freedom for the prior distribution are estimated from the data.

The effect of this function is to smooth or shrink the variances towards a common value.
The smoothed variances have a smaller expected mean square error to the true variances than do the sample variances themselves.

This function is called by \code{eBayes}, but beware a possible confusion with the output from that function.
The values \code{var.prior} and \var{var.post} output by \code{squeezeVar} correspond to the quantities \code{s2.prior} and \code{s2.post} output by \code{eBayes}, whereas \code{var.prior} output by \code{eBayes} relates to a different parameter.
}

\value{
A list with components
\item{var.post}{numeric vector of posterior variances.}
\item{var.prior}{location of prior distribution. A vector if \code{covariate} is non-\code{NULL}, otherwise a scalar.}
\item{df.prior}{degrees of freedom of prior distribution.}
}

\author{Gordon Smyth}

\references{
Smyth, G. K. (2004). Linear models and empirical Bayes methods for
     assessing differential expression in microarray experiments.
     \emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3},
     No. 1, Article 3. \url{http://www.bepress.com/sagmb/vol3/iss1/art3}
}

\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\examples{
s2 <- rchisq(20,df=5)/5
squeezeVar(s2, df=5)
}

\keyword{htest}
