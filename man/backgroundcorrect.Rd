\name{backgroundCorrect}
\alias{backgroundCorrect}
\title{Correct Intensities for Background}
\description{
Apply background correction to microarray expression intensities.
}
\usage{
backgroundCorrect(RG, method="subtract")
}
\arguments{
  \item{RG}{an \code{\link[limma:rglist]{RGList}} object or a unclassed list containing the same components as an \code{RGList}}
  \item{method}{character string specifying correction method.  Possible values are \code{"none"}, \code{"subtract"}, \code{"half"} or \code{"minimum"}.}
}
\details{
If \code{method="none"} then the corrected intensities are equal to the foreground intensities, i.e., the background intensities are treated as zero.
If \code{method="subtract"} then this function simply subtracts the background intensities from the foreground intensities which is the usual background correction method.
If \code{method="half"} then any intensity which is less than 0.5 after background subtraction is reset to be equal to 0.5.
If \code{method="minimum"} then any intensity which is zero or negative after background subtraction is set equal to half the minimum of the positive corrected intensities for that array.

Background correction (background subtraction) is also performed by the \code{\link{normalizeWithinArrays}} method for \code{RGList} objects, so it is not necessary to call \code{backgroundCorrect} directly unless one wants to use a method other than simple subtraction.
Calling \code{backgroundCorrect} before \code{normalizeWithinArrays} will over-ride the default background correction.
}
\value{
An \code{RGList} object in which components \code{R} and \code{G} are background corrected
and components \code{Rb} and \code{Gb} are removed.
}
\references{
Yang, Y. H., Buckley, M. J., Dudoit, S., and Speed, T. P. (2002). Comparison of methods for image analysis on cDNA microarray data. Journal of Computational and Graphical Statistics 11, 108-136.

Yang, Y. H., Buckley, M. J., and Speed, T. P. (2001). Analysis of microarray images. Briefings in Bioinformatics 2, 341-349.
}
\author{Gordon Smyth}
\examples{
RG <- new("RGList", list(R=c(1,2,3,4),G=c(1,2,3,4),Rb=c(2,2,2,2),Gb=c(2,2,2,2)))
backgroundCorrect(RG)
backgroundCorrect(RG, method="half")
backgroundCorrect(RG, method="minimum")
}
\seealso{
	\code{\link{read.maimages}}, \code{\link{normalizeWithinArrays}}
}
\keyword{models}
