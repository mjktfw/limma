\name{nec}
\alias{nec}
\title{NormExp by Control background correction (NEC)}
\description{Perform normexp background correction using negative control probes.}
\usage{
nec(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE,detection.p="Detection")
}
\arguments{
  \item{x}{object of class \code{\link{EListRaw-class}} or \code{matrix} containing raw intensities for regular and control probes from a series of microarrays.}
  \item{status}{character vector giving probe types.}
  \item{negctrl}{character string identifier for negative control probes.}
  \item{regular}{character string identifier for regular probes.}
  \item{offset}{numeric value added to the intensities after background correction.}
  \item{robust}{logical. Should robust estimators be used for the background mean and standard deviation?}
  \item{detection.p}{a character string giving the name of the component which contains detection p value information in \code{x} or a numeric matrix giving detection p values, \code{Detection} by default}
  }
\details{
When control data are available, this function calls \code{\link{normexp.fit.control}} to estimate the parameters required by normal+exponential(normexp) convolution model with the help of negative control probes, followed by \code{\link{normexp.signal}} to perform the background correction.
If \code{x} contains background intensities \code{x$Eb}, then these are first subtracted from the foreground intensities, prior to normexp background correction.
After background correction, an \code{offset} is added to the data.

When control data are not available, this function calls \code{\link{normexp.fit.detection.p}} to estimate the normexp parameters.
\code{\link{normexp.fit.detection.p}} infers negative control probe intensities from regular probes by taking advantage of their detection p value information.

For more descriptions to parameters \code{x}, \code{status}, \code{negctrl}, \code{regular} and \code{detection.p}, please refer to functions \code{\link{normexp.fit.control}}, \code{\link{normexp.fit.detection.p}} and \code{\link{read.ilmn}}.
}
\value{
An \code{\link{EListRaw-class}} or matrix object containing background-corrected intensities. Intensities are at raw scale. Control probes are kept. 
}

\references{
Shi W, Oshlack A and Smyth GK (2010).
Optimizing the noise versus bias trade-off for Illumina Whole Genome Expression BeadChips.
\emph{Nucleic Acids Research} 38, e204.
\url{http://nar.oxfordjournals.org/content/38/22/e204}
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
  An overview of background correction functions is given in \link{04.Background}.
  
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

  \code{\link{normexp.fit.control}} estimates the parameters in the normal+exponential convolution model using the negative control probes.

  \code{\link{normexp.fit.detection.p}} estimates the parameters in the normal+exponential convolution model using negative control probe intensities inferred from regular probes by using their detection p values information.

  
  \code{\link{normexp.fit}} estimates parameters in the normal+exponential convolution model using a saddle-point approximation or other methods.

  \code{\link{neqc}} performs normexp background correction and quantile normalization aided by control probes.  
}

\examples{
\dontrun{
x <- read.ilmn(files="sample probe profile.txt",ctrlfiles="control probe profile.txt")
x.bg <- nec(x)
}
}

\keyword{models}
