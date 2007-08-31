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
  \item{par}{numeric vector containing the parameters of the Normal+Exp distribution, see \code{normexp.fit}}
  \item{x}{numeric vector of (background corrected) intensities}
}
\details{
In general the vector \code{normmean} is computed conditional on background at each spot.
}
\value{
Numeric vector containing adjusted intensities.
}

\references{
McGee, M., and Chen, Z. (2006).
Parameter estimation for the exponential-normal convolution model for background correction of Affymetrix GeneChip data.
\emph{Stat Appl Genet Mol Biol}, 5(1), Article 24.

Ritchie, M. E., Silver, J., Oshlack, A., Silver, J., Holmes, M., Diyagama, D., Holloway, A., and Smyth, G. K. (2007).
A comparison of background correction methods for two-colour microarrays.
\emph{Bioinformatics}
\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btm412}
}

\author{Gordon Smyth}
\seealso{
\code{\link{normexp.fit}}

An overview of background correction functions is given in \code{\link{04.Background}}.
}
\examples{# See normexp.fit}
\keyword{models}
