\title{MA-Plot}
\name{plotMA}
\alias{plotMA}
\description{
Creates an MA-plot with color coding for various sort of control spots.
}
\usage{
plotMA(MA,array=1,pch=16,status=NULL,
       values=c("gene","blank","buffer","utility","negative","calibration","ratio"),
       col=c("black","yellow","orange","pink","brown","blue","red"),
       cex=c(0.1,0.6,0.6,0.6,0.6,0.6,0.6))
}
\arguments{
  \item{MA}{list with components \code{M} containing log-ratios and \code{A} containing average intensities}
  \item{array}{integer giving the array to be plotted. Corresponds to columns of \code{M} and \code{A}.}
  \item{pch}{plotting character}
  \item{status}{character vector giving the control status of each spot on the array}
  \item{values}{character vector giving unique values of \code{status} corresponding to control states of interest}
  \item{col}{vector of colors, of the same length as \code{values}}
  \item{cex}{numeric vector of the same length as \code{values} giving sizes for plot symbols}
}
\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{6.Diagnostics}.
}
\keyword{hplot}
