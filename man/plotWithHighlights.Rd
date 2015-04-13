\name{plotWithHighlights}
\alias{plotWithHighlights}
\title{Scatterplot With Highlighting of Special Points}
\description{
Creates scatterplot, with optional size and color coding for points of special interest.
This is the engine for \code{plotMD} and \code{plotMA}.
}
\usage{
plotWithHighlights(x, y, status = NULL, values = NULL,
                   hi.pch = 16, hi.col = NULL, hi.cex = 1, legend = "topleft",
                   bg.pch = 16, bg.col = "black", bg.cex = 0.3, \dots)
}
\arguments{
  \item{x}{numeric vector.}

  \item{y}{numeric vector.}

  \item{status}{character vector giving the control status of each point, of same length as \code{x} and \code{y}.
  If \code{NULL}, then all points are plotted in the default color, symbol and size.}

  \item{values}{character vector giving values of \code{status} to be highlighted on the plot.
  Defaults to unique values of \code{status} in decreasing order of frequency, with the most frequent value set as the background value.
  Ignored if there is no \code{status} vector.}

  \item{hi.pch}{vector of plotting characters for highlighted points, either of unit length or of same length as \code{values}.
  Ignored is there is no \code{status} vector.}

  \item{hi.col}{vector of colors for highlighted points, either of unit length or of same length as \code{values}.
  Defaults to \code{1+1:length(values)}.
  Ignored if there is no \code{status} vector.}

  \item{hi.cex}{numeric vector of plot symbol expansions for highlighted points, either of unit length or of same length as \code{values}.
  Ignored if there is no \code{status} vector.}

  \item{legend}{character string giving position to place legend.
  See \code{\link{legend}} for possible values.
  Can also be logical, with \code{FALSE} meaning no legend.
  Ignored if there is no \code{status} vector.}

  \item{bg.pch}{plotting character for background (non-highlighted) points.}

  \item{bg.col}{color for background (non-highlighted) points.}

  \item{bg.cex}{plot symbol expansion for background (non-highlighted) points.}

  \item{\dots}{other arguments are passed to \code{plot}.}
}

\details{
This function produces a scatterplot in which the highlighted points are, by default, larger and colored compared to background points.

The \code{status} vector might be used specify the control status of each point, for example \code{"gene"}, \code{"ratio control"}, \code{"house keeping gene"}, \code{"buffer"} and so on.
The vector can be computed using the function \code{\link{controlStatus}} and a spot-types file.
Alternatively, the \code{status} vector could indicate statistically significant points.
However the function may be used to highlight any subset of spots.

By default, the most common value of \code{status} is assumed to correspond to a background point.

See \code{\link[graphics]{points}} for possible values for \code{pch}, \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43,
doi: 10.1093/nar/gkv007.
}

\author{Gordon Smyth}

\examples{
x <- runif(1000,4,16)
status <- rep(c(0,-1,1),c(950,40,10))
y <- status + rnorm(1000,sd=0.2)
plotWithHighlights(x,y,status=status)
}

\seealso{
\code{\link{plotMD}}

An overview of diagnostic plots available in LIMMA is given in \link{09.Diagnostics}.
}
