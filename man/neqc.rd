\name{neqc}
\alias{neqc}
\title{NormExp and Quantile by Control (NEQC)}
\description{Perform normexp background correction and quantile normalization using control probes.}
\usage{
neqc(x, status=NULL, negctrl="negative", regular="regular", offset=16, robust=FALSE, ...)
}
\arguments{
  \item{x}{object of class \code{\link{EListRaw-class}} or \code{matrix} containing raw intensities for regular and control probes from a series of microarrays.}
  \item{status}{character vector giving probe types.}
  \item{negctrl}{character string identifier for negative control probes.}
  \item{regular}{character string identifier for regular probes.}
  \item{offset}{numeric value added to the intensities after background correction.}
  \item{robust}{logical. Should robust estimators be used for the background mean and standard deviation?}
  \item{...}{any other arguments are passed to \code{normalizeBetweenArrays.}}
  }
\details{
This function calls \code{\link{nec}} to carry out normexp background correction
aided by negative control probes.
\code{\link{nec}} will add an offset to the data after background correction.
If \code{x} contains background intensities \code{x$Eb}, then these are first subtracted from the foreground intensities, prior to normexp background correction.

The background-corrected intensities will then be quantile normalized, including control probes.
Finally the intensities are log2 transformed and the control probes are removed.

For more descriptions to parameters \code{x}, \code{status}, \code{negctrl} and \code{regular}, please refer to functions \code{\link{normexp.fit.control}} and \code{\link{read.ilmn}}.
}
\value{
An \code{\link{EList-class}} or matrix object containing normalized log2 expression values. Control probes are removed.
}

\references{
Shi W, Oshlack A and Smyth GK (2010). Optimizing the noise versus bias trade-off for Illumina Whole Genome Expression BeadChips. Nucleic Acids Research, October 6. [Epub ahead of print]. PMID: 20929874 
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
  
  An overview of background correction functions is given in \link{04.Background}.

  \code{\link{nec}} performs normexp by control background correction. 
  
  \code{\link{normexp.fit.control}} estimates the parameters in the normal+exponential convolution model using the negative control probes.
  
  \code{\link{normexp.fit}} estimates parameters in the normal+exponential convolution model using a saddle-point approximation or other methods.
}

\examples{
\dontrun{
x <- read.ilmn(files="sample probe profile.txt",ctrlfiles="control probe profile.txt")
y <- neqc(x)
}
}

\keyword{models}
