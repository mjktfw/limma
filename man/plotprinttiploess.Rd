\title{MA Plots by Print-Tip Group}
\name{plotPrintTipLoess}
\alias{plotPrintTipLoess}
\description{
Creates a coplot giving MA-plots with loess curves by print-tip groups.
}
\usage{
plotPrintTipLoess(MA,layout,array=1)
}
\arguments{
  \item{MA}{list with components \code{M} containing log-ratios and \code{A} containing average intensities}
  \item{layout}{a list specifying the number of tip rows and columns and the number of spot rows and columns printed by each tip}
  \item{array}{integer giving the array to be plotted. Corresponds to columns of \code{M} and \code{A}.}
}
\value{A plot is created on the current graphics device}
\author{Gordon Smyth}
\keyword{hplot}
