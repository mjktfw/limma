\name{normexp.signal}
\alias{normexp.signal}
\title{Expected Signal Given Observed Foreground Under Normal+Exp Model}
\description{
Adjust foreground intensities for observed background using Normal+Exp Model.
This function is called by \code{backgroundCorrect} and is not normally called directly by the user.
}
\usage{
normexp.signal(par, x)
}
\arguments{
  \item{par}{numeric vector containing the parameters of the Normal+Exp distribution, see \code{\link{normexp.fit}} for details.}
  \item{x}{numeric vector of (background corrected) intensities}
}
\details{
In general the vector \code{normmean} is computed conditional on background at each spot.
}
\value{
Numeric vector containing adjusted intensities.
}

\references{
Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
\emph{Bioinformatics}
\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btm412}

Silver, JD, Ritchie, ME, and Smyth, GK (2009).
Microarray background correction: maximum likelihood estimation for the normal-exponential convolution.
\emph{Biostatistics} 10, 352-363.
\url{http://biostatistics.oxfordjournals.org/cgi/content/abstract/kxn042}
}

\author{Gordon Smyth}
\seealso{
\code{\link{normexp.fit}}

An overview of background correction functions is given in \code{\link{04.Background}}.
}
\examples{# See normexp.fit}
\keyword{models}
