\name{normexp.fit}
\alias{normexp.fit}
\alias{normexp.fit.C}
\title{Fit Normal+Exp Convolution Model to Observed Intensities}
\description{
Fit normal+exponential convolution model to observed intensities.
The normal part represents the background and the exponential represents the signal intensities.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.fit(x, method="saddle", n.pts=NULL, trace=FALSE)
normexp.fit.C(x, method="saddle", n.pts=NULL, trace=FALSE)
}
\arguments{
  \item{x}{numeric vector of (background corrected) intensities}
  \item{method}{method used to estimate the three parameters. Choices for \code{normexp.fit} are \code{"saddle"}, \code{"neldermead"}, \code{"bfgs"}, \code{"rma"} and \code{"mcgee"}.
  Additional choices for \code{normexp.fit.C} are \code{"nlminb"} and \code{"nlminblog"}.}
  \item{n.pts}{number of quantiles of \code{x} to use for the fit. If \code{NA} then all values of \code{x} will be used.}
  \item{trace}{logical, if \code{TRUE}, tracing information on the progress of the optimization is given.}
}
\details{
\code{normexp.fit} is a pure R function while \code{normexp.fit.C} call C source code.
The two functions should give identical values but the latter is several times as fast.

This function uses maximum likelihood estimation to fit a model to the background-corrected intensities.
The model is a $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the background corrected intensities.

This is essentially the same model which is used by \code{bg.correct.rma} in the affy package.
The difference is that the parameters are estimated by maximum likelihood and that the estimated background is subtracted before applying the model-based background.

For speed and numerical stability, the likelihood function is approximated by default using a saddle-point approximation.
Other choices are \code{"neldermead"} or \code{"bfgs"} which attempt to maximize the exact log-likelihood, which may be numerically unstable.
\code{"rma"} uses the kernel estimation method from the RMA algorithm of affy package.
\code{"mcgee"} uses the RMA-75 method from McGee and Chen (2006).
\code{"nlminb"} and \code{"nlminblog"} compute the maximum likelihood estimates.

If the number of rows of \code{x} is large, computation time can be saved by setting \code{n.pts} to a value less than the total number of probes, for example \code{n.pts=2^14}.
}
\value{
A list containing the components
  \item{par}{numeric vector giving estimated values of $\mu$, $\log(\sigma)$ and $\log\alpha$}
  \item{m2loglik}{numeric scalar giving minus twice the log-likelihood}
  \item{convergence}{integer code indicating successful convergence or otherwise of the optimization. See \code{optim}.}
}
\author{Gordon Smyth and Jeremy Silver}

\references{
McGee, M., and Chen, Z. (2006).
Parameter estimation for the exponential-normal convolution model for background correction of Affymetrix GeneChip data.
\emph{Stat Appl Genet Mol Biol}, 5(1), Article 24.

Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
\emph{Bioinformatics}
\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btm412}
}

\seealso{
\code{\link{normexp.signal}}

An overview of background correction functions is given in \code{\link{04.Background}}.
}
\examples{
f <- c(2,3,1,10,3,20,5,6)
b <- c(2,2,2,2,2,2,2,2)
out <- normexp.fit(f-b)
normexp.signal(out$par, x=f-b)
}
\keyword{models}
