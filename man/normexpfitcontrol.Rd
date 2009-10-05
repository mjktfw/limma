\name{normexp.fit.control}
\alias{normexp.fit.control}
\title{Normexp Model Parameter Estimation Aided by Negative Controls}
\description{The mean and log-standard-deviation of the background-normal part of the normexp+exponential convolution model is estimated as the mean and log-standard deviation of intensities from negative control probes. The log-mean of the signal-exponential part is estimated as the log of the difference between signal mean and background mean.}
\usage{
normexp.fit.control(x, status=NULL, negctrl="negative")
}
\arguments{
  \item{x}{ object of class \code{EListRaw-class} or \code{matrix} containing raw intensities for regular and control probes for a series of microarrays}
  \item{status}{ character vector giving probe types.}
  \item{negctrl}{ character string giving the type identifier of negative control probes.}
  }
\details{
\code{x} has to contain raw expression intensities from both regular probes and negative control probes.

The probe type information for an object of \code{\link{EListRaw-class}} is normally saved in the \code{Status} column of its \code{genes} component.
However, it will be overriden by the \code{status} parameter if it is explicitly provided to this function.
If \code{x} is a \code{matrix} object, the probe type information has to be provided through the \code{status} parameter of this function.
Regular probes have the status \code{regular}.
Negative control probes have the status indicated by \code{negctrl}, which is \code{negative} by default.

This function estimates parameters of the normal+exponential convolution model with the help of negative control probes.
The mean and log-standard-deviation of the background-normal part of the normexp+exponential convolution model are estimated as the mean and log-standard deviation of intensities from negative control probes respectively.
The log-mean of the signal-exponential part is estimated as the log of the difference between signal mean and background mean.
The signal mean is simply the mean of intensities from regular probes.
}
\value{
A matrix containing estimated parameters with rows being arrays and with columns being parameters.
Column names are \code{mu}, \code{logsigma} and \code{logalpha}.
}

\references{
Wei Shi and Gordon K Smyth. Normalizing Illumina Whole Genome Expression BeadChips. In preparation.
}

\author{Wei Shi and Gordon Smyth}

\seealso{
\code{\link{neqc}} calls this function to get the parameters of the normal+exponential convolution model and then calls \code{\link{normexp.signal}} to perform the background correction.

\code{\link{normexp.fit}} estimates parameters in the normal+exponential convolution model using a saddle-point approximation or other mothods.

An overview of background correction functions is given in \code{\link{04.Background}}.
}

\examples{
\dontrun{
x <- read.profile(files="sample probe profile", ctrlfiles="control probe profile")
par <- normexp.fit.control(x)
}
}

\keyword{models}