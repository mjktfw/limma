\name{backgroundCorrect}
\alias{backgroundCorrect}
\title{Correct Intensities for Background}
\description{
Apply background correction to microarray expression intensities.
}
\usage{
backgroundCorrect(RG)
}
\arguments{
  \item{RG}{an \code{\link[limma:rglist]{RGList}} object or a unclassed list containing the same components as an \code{RGList}}
}
\details{
At present this function simply subtracts the background intensities from the foreground intensities.
Other correction methods are planned.

Background correction is also performed by the \code{\link[limma:normalizeWithinArrays-methods]{normalizeWithinArrays}} method for \code{RGList} objects, so at this time it is not necessary to call \code{backgroundCorrect} directly.
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
\seealso{
	\code{\link{read.maimages}}, \code{\link{normalizeWithinArrays}}
}
\keyword{models}
