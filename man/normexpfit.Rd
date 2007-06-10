\name{normexp.fit}
\alias{normexp.fit}
\title{Fit Normal+Exp Convolution Model to Observed Intensities}
\description{
Fit normal+exponential convolution model to observed intensities.
The normal part represents the background and the exponential represents the signal intensities.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.fit(x, n.pts=2^10, trace=FALSE)
}
\arguments{
  \item{x}{numeric vector of (background corrected) intensities}
  \item{n.pts}{number of quantiles of \code{x} to use for the fit. If \code{NA} then all values of \code{x} will be used.}
  \item{trace}{logical, if \code{TRUE}, tracing information on the progress of the optimization is given.}
}
\details{
This function uses maximum likelihood estimation to fit a model to the foreground and background intensities.
The model is a $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the background corrected intensities.

This is essentially the same model which is used by \code{bg.correct.rma} in the affy package.
The difference is that the parameters are estimated by maximum likelihood and that the estimated background is subtracted before applying the model-based background.

For speed and numerical stability, the likelihood function is approximated using a saddle-point approximation.

By default, about 1000 quantiles of \code{x} are used of the fit.
This limits computation time, and also allows a few outliers to be ignored at the top and bottom of the distribution if the length of \code{n.pts} is less than \code{length(x)}.
}
\value{
A list containing the components
  \item{par}{numeric vector giving estimated values of $\mu$, $\log(\sigma)$ and $\log\alpha$}
  \item{m2loglik}{numeric scalar giving minus twice the log-likelihood}
  \item{convergence}{integer code indicating successful convergence or otherwise of the optimization. See \code{optim}.}
}
\author{Gordon Smyth and Jeremy Silver}

\references{
Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
Submitted.
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
