\title{MA Plots by Print-Tip Group}
\name{plotPrintTipLoess}
\alias{plotPrintTipLoess}
\description{
Creates a coplot giving MA-plots with lowess curves by print-tip groups.
}
\usage{
plotPrintTipLoess(MA,layout,array=1,span=0.3,...)
}
\arguments{
  \item{MA}{list with components \code{M} containing log-ratios and \code{A} containing average intensities}
  \item{layout}{a list specifying the number of tip rows and columns and the number of spot rows and columns printed by each tip}
  \item{array}{integer giving the array to be plotted. Corresponds to columns of \code{M} and \code{A}.}
  \item{span}{span of window for \code{lowess} curve}
  \item{...}{other arguments passed to \code{panel.smooth}}
}
\value{A plot is created on the current graphics device}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{6.Diagnostics}.
}
\author{Gordon Smyth}
\keyword{hplot}
