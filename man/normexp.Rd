\name{normexp}
\alias{normexp.m2loglik}
\alias{normexp.m2loglik.saddle}
\alias{normexp.grad}
\title{Normal + Exponential Log-Likelihood}
\description{
Marginal log-likelihood of foreground values for the normal + exponential convolution model and its derivatives.
These functions are called by \code{normexp.fit} and are not normally called directly by the user.
}
\usage{
normexp.m2loglik(par,x)
normexp.m2loglik.saddle(par,x)
normexp.grad(par,x)
}
\arguments{
  \item{par}{numeric vector of parameters}
  \item{x}{numeric vector of (background corrected) intensities}
}
\details{
The parameter vector \code{par} holds the normal mean, the normal log-standard deviation and the exponential mean.

\code{normexp.m2loglik} computes minus twice the log-likelihood, and \code{normexp.grad} it is derivative, based on the $normal(\mu,\sigma^2)+exponential(\alpha)$ convolution model for the intensities.
The elements of \code{par} are $\mu$, $\log(\sigma)$ and $\log(\alpha)$.

\code{normexp.m2loglik} is the saddle-point approximation to the log-logelihood, which is generally prefered because it is numerically more stable.
}
\value{
\code{normexp.m2loglik} returns a numeric scalar holding minus-twice the log-likelihood.
\code{normexp.grad} returns a numeric vector holding the derivatives with respect to the elements of \code{par}.
}
\author{Jeremy Silver and Gordon Smyth}
\seealso{
An overview of background correction functions is given in \code{\link{04.Background}}.
}
\keyword{models}
