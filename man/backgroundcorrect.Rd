\name{backgroundCorrect}
\alias{backgroundCorrect}
\title{Correct Intensities for Background}
\description{
Background correct microarray expression intensities.
}
\usage{
backgroundCorrect(RG, method="subtract", printer=RG$printer)
}
\arguments{
  \item{RG}{an \code{\link[limma:rglist]{RGList}} object or a unclassed list containing the same components as an \code{RGList}}
  \item{method}{character string specifying correction method.  Possible values are \code{"none"}, \code{"subtract"}, \code{"half"}, \code{"minimum"}, \code{"movingmin"} or \code{"edwards"}.}
\item{printer}{a list containing printer layout information, see \code{\link{PrintLayout-class}}}
}
\details{
If \code{method="none"} then the corrected intensities are equal to the foreground intensities, i.e., the background intensities are treated as zero.
If \code{method="subtract"} then this function simply subtracts the background intensities from the foreground intensities which is the usual background correction method.
If \code{method="half"} then any intensity which is less than 0.5 after background subtraction is reset to be equal to 0.5.
If \code{method="minimum"} then any intensity which is zero or negative after background subtraction is set equal to half the minimum of the positive corrected intensities for that array.
If \code{method="movingmin"} then the background estimates are replaced with the minimums of the backgrounds of the spot and its eight neighbors, i.e., the background is replaced by a moving minimum of 3x3 grids of spots.
If \code{method="edwards"} a log-linear interpolation method is used to adjust lower intensities as in Edwards (2003).

Background correction (background subtraction) is also performed by the \code{\link{normalizeWithinArrays}} method for \code{RGList} objects, so it is not necessary to call \code{backgroundCorrect} directly unless one wants to use a method other than simple subtraction.
Calling \code{backgroundCorrect} before \code{normalizeWithinArrays} will over-ride the default background correction.
}
\value{
An \code{RGList} object in which components \code{R} and \code{G} are background corrected
and components \code{Rb} and \code{Gb} are removed.
}
\references{
Edwards, D. E. (2003). Non-linear normalization and background correction in one-channel cDNA microarray studies
\emph{Bioinformatics} 19, 825-833. 

Yang, Y. H., Buckley, M. J., Dudoit, S., and Speed, T. P. (2002). Comparison of methods for image analysis on cDNA microarray data. \emph{Journal of Computational and Graphical Statistics} 11, 108-136.

Yang, Y. H., Buckley, M. J., and Speed, T. P. (2001). Analysis of microarray images. \emph{Briefings in Bioinformatics} 2, 341-349.
}
\author{Gordon Smyth}
\examples{
RG <- new("RGList", list(R=c(1,2,3,4),G=c(1,2,3,4),Rb=c(2,2,2,2),Gb=c(2,2,2,2)))
backgroundCorrect(RG)
backgroundCorrect(RG, method="half")
backgroundCorrect(RG, method="minimum")
}
\seealso{
An overview of normalization and background correction functions is given in \code{\link{4.Normalization}}.
}
\keyword{models}
