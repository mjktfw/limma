\title{Kooperberg Model-Based Background Correction}
\name{kooperberg}
\alias{kooperberg}
\description{
This function uses a Bayesian model to background correct
data from a series of microarray experiments.
It is currently implemented only for GenePix data.
}

\usage{
kooperberg(names, fg="mean", bg="median")
}
\arguments{
\item{names}{character vector giving the names of data.frames containing GenePix data}
\item{fg}{character string giving foreground estimator.
Choices are \code{"mean"} or \code{"median"}.}
\item{bg}{character string giving foreground estimator.
Choices are \code{"mean"} or \code{"median"}.}
}

\details{
\code{bayesianAdjustedRG} serially extracts the foreground and background intensities, standard
deviations and number of pixels from GenePix data frames.

This information provides empirical estimates of the model parameters
described in equation 2 of Kooperberg et al (2003).
Note that the 'a' parameters in equations 3 and 4 have been omitted for simplicity.

Choices are \code{"mean"} for \code{F635.Mean} (red) and
\code{F532.Mean} (green) or \code{"median"} for \code{F635.Median}
and \code{F532.Median}.}
\item{bg}{character string giving foreground estimator.
Choices are \code{"mean"} for \code{B635.Mean} (red) and
\code{B532.Mean} (green) or \code{"median"} for \code{B635.Median}
and \code{B532.Median}.}

}

\value{
A list containing the components
\item{R}{matrix containing the background adjusted intensities for
the red channel for each spot for each array}
\item{G}{matrix containing the background adjusted intensities for the green channel for each spot for each array}
}

\author{Matthew Ritchie}

\references{
Kooperberg, C., Fazzio, T. G., Delrow, J. J., and Tsukiyama, T. (2002)
Improved background correction for spotted DNA microarrays.
\emph{Journal of Computational Biology} \bold{9}, 55-66.
}
	
\seealso{
\link{4.Normalization} gives an overview of normalization and background correction functions defined in the LIMMA package.
}

\examples{
\dontrun{genepixFiles <- dir(pattern="\\\\.gpr") # get the names of the GenePix image analysis output files in the current directory
read.series(genepixFiles, suffix=NULL, skip=26, sep="\t") # read in GenePix files
RGmodel <- bayesianAdjustedRG(genepixFiles) # perform model based background correction for GenePix data}
}

\keyword{models}
