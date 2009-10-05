\name{neqc}
\alias{neqc}
\title{NormExp and Quantile by Control (NEQC)}
\description{Perform normexp background correction and quantile normalization using control probes.}
\usage{
neqc(x, status=NULL, negctrl="negative", regular="regular", offset=50, ...)
}
\arguments{
  \item{x}{ object of class \code{\link{EListRaw-class}} or \code{matrix} containing raw intensities for regular and control probes from a series of microarrays.}
  \item{status}{ character vector giving probe types.}
  \item{negctrl}{ character string giving the \code{status} identifier of negative control probes. Default is \code{negative}.}
  \item{regular}{ character string giving the \code{status} identifier of regular probes. Default is \code{regular}.}
  \item{offset}{ numeric value added to the intensities after background correction.}
  \item{...}{ any other parameters are passed to \code{normalizeBetweenArrays.}}
  }
\details{
This function calls the function \code{\link{normexp.fit.control}} to estimate the parameters required by normal+exponential convolution model with the help of negative control probes.
\code{\link{normexp.signal}} is then called to background correct the raw data.
An \code{offset} is added to the data after the background correction.
This function will then call the function \code{\link{normalizeBetweenArrays}} to perform quantile between-array normalization and log2 transformation.

For more descriptions to parameters \code{x}, \code{status}, \code{negctrl} and \code{regular}, please refer to functions \code{\link{normexp.fit.control}} and \code{\link{read.ilmn}}.
}
\value{
An \code{\link{EList-class}} object containing normalized log2 expression values. Control probes are removed.
}

\references{
Wei Shi and Gordon K Smyth. Normalizing Illumina Whole Genome Expression BeadChips. In preparation.
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
  An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
  
  An overview of background correction functions is given in \link{04.Background}.
  
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