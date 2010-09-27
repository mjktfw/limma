\name{nec}
\alias{nec}
\title{NormExp by Control background correction (NEC)}
\description{Perform normexp background correction using negative control probes.}
\usage{
nec(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE)
}
\arguments{
  \item{x}{object of class \code{\link{EListRaw-class}} or \code{matrix} containing raw intensities for regular and control probes from a series of microarrays.}
  \item{status}{character vector giving probe types.}
  \item{negctrl}{character string identifier for negative control probes.}
  \item{regular}{character string identifier for regular probes.}
  \item{offset}{numeric value added to the intensities after background correction.}
  \item{robust}{logical. Should robust estimators be used for the background mean and standard deviation?}
  }
\details{
This function calls \code{\link{normexp.fit.control}} to estimate the parameters required by normal+exponential convolution model with the help of negative control probes, followed by \code{\link{normexp.signal}} to perform the background correction.
If \code{x} contains background intensities \code{x$Eb}, then these are first subtracted from the foreground intensities, prior to normexp background correction.
After background correction, an \code{offset} is added to the data.

For more descriptions to parameters \code{x}, \code{status}, \code{negctrl} and \code{regular}, please refer to functions \code{\link{normexp.fit.control}} and \code{\link{read.ilmn}}.
}
\value{
An \code{\link{EListRaw-class}} or matrix object containing background-corrected intensities. Intensities are at raw scale. Control probes are kept. 
}

\references{
Shi W, Oshlack A and Smyth GK (2010). Optimizing the noise versus bias trade-off for Illumina Whole Genome Expression BeadChips. Nucleic Acids Research, In Press 
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
  An overview of background correction functions is given in \link{04.Background}.
  
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.

  \code{\link{normexp.fit.control}} estimates the parameters in the normal+exponential convolution model using the negative control probes.
  
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
