\title{Write MA-Plots to Files}
\name{plotMA3by2}
\alias{plotMA3by2}
\description{
Write MA-plots to files in PNG format, six plots to a file in a 3 by 2 grid arrangement.
}
\usage{
plotMA3by2(MA, prefix="MA", path=NULL, main=colnames(MA), zero.weights=FALSE, common.lim=TRUE, device="png", ...)
}
\arguments{
  \item{MA}{an \code{MAList} or \code{RGList} object, or any list with components \code{M} containing log-ratios and \code{A} containing average intensities}
  \item{prefix}{character string giving prefix to attach to file names}
  \item{path}{character string specifying directory for output files}
  \item{main}{character vector giving titles for plots}
  \item{zero.weights}{logical, should points with non-positive weights be plotted}
  \item{common.lim}{logical, should all plots on a page use the same axis limits}
  \item{device}{device driver for the plot. Choices are \code{"png"}, \code{"jpeg"}, \code{"pdf"}, \code{"postscript"}.}
  \item{...}{any other arguments are passed to \code{plotMA}}
}

\details{
This function writes a series of graphic files to disk.
Each file contains six MA-plots in three rows and two columns.
The layout is optimized for A4-sized paper.

The graph format can be \code{"png"} or \code{"jpeg"}, which are screen-resolution formats, or \code{"pdf"} or \code{"postscript"}, which are loss-less formats.
\code{"png"} is not available on every R platform.
Note that \code{"pdf"} or \code{"postscript"} may produce very large files.
}

\value{
No value is returned, but one or more files are written to the working directory.
The number of files is determined by the number of columns of \code{MA}.
}

\author{Gordon Smyth}

\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}
