\name{backgroundCorrect}
\alias{backgroundCorrect}
\title{Background Correct}
\description{
Apply background correction to microarray expression intensities.
}
\usage{
backgroundCorrect(RG)
}
\arguments{
  \item{RG}{A list containing the components \code{R}, \code{G}, \code{Rb} and \code{Gb} containing the foreground and background intensities for each spot on each array.}
}
\details{
At present this function simply subtracts the background intensities from the foreground intensities.
Other correction methods are planned.

Background correction is also performed by the \code{\link{normalizeWithinArrays.RGList}} normalize method, so it is not necessary to call \code{backgroundCorrect} explicitly.
}
\value{
A list the same as \code{RG} except that components \code{R} and \code{G} are background corrected
and components \code{Rb} and \code{Gb} are removed.
}
\references{
Yang, Y. H., Buckley, M. J., Dudoit, S., and Speed, T. P. (2002). Comparison of methods for image analysis on cDNA microarray data. Journal of Computational and Graphical Statistics 11, 108-136.

Yang, Y. H., Buckley, M. J., and Speed, T. P. (2001). Analysis of microarray images. Briefings in Bioinformatics 2, 341-349.
}
\author{Gordon Smyth}
\seealso{
	\code{\link{read.maimages}}, \code{\link{normalizeWithinArrays.RGList}}
}
\keyword{models}
