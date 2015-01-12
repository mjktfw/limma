\title{plotWithHighlights}
\name{plot.MArrayLM}
\alias{plot.RGList}
\alias{plot.MAList}
\alias{plot.EList}
\alias{plot.MArrayLM}
\alias{plotWithHighlights}
\description{
Plot data or model fit objects.
These represent the object by creating an MA-plot, with optional color coding for control spots.
}
\usage{
plotWithHighlights(x, y, status = NULL, values = NULL, pch = 16, col = NULL, cex = 1,
     legend = "topleft", pch.bg = 16, col.bg = "black", cex.bg = 0.3, \dots)
\method{plot}{EList}(x, y, array = 1, xlab = "Average log-expression",
     ylab = "Expression log-ratio (this sample vs others)", main = colnames(x)[array],
     status=x$genes$Status, zero.weights = FALSE, \dots)
\method{plot}{RGList}(x, y, array = 1, xlab = "A", ylab = "M", main = colnames(x)[array],
     status=x$genes$Status, zero.weights = FALSE, \dots)
\method{plot}{MAList}(x, y, array = 1, xlab = "A", ylab = "M", main = colnames(x)[array],
     status=x$genes$Status, zero.weights = FALSE, \dots)
\method{plot}{MArrayLM}(x, y, coef = ncol(x), xlab = "Average log-expression", ylab = "log-fold-change",
     main = colnames(x)[coef], status = x$genes$Status, zero.weights = FALSE, \dots)
}
\arguments{
  \item{x}{an \code{RGList}, \code{MAList}, \code{EList} or \code{MArrayLM} object.}
  \item{y}{not used.}
  \item{array}{integer giving the array to be plotted (if \code{x} is an \code{RGList}, \code{MAList} or \code{EList} object).}
  \item{coef}{integer giving the linear model coefficient to be plotted.}
  \item{xlab}{character string giving label for x-axis}
  \item{ylab}{character string giving label for y-axis}
  \item{main}{character string giving title for plot}
  \item{status}{character vector giving the control status of each spot on the array, of same length as the number of rows of \code{MA$M}.
  If \code{NULL}, then all points are plotted in the default color, symbol and size.}
  \item{values}{character vector giving values of \code{status} to be highlighted on the plot.
  Defaults to unique values of \code{status} in decreasing order of frequency, with the most frequent value set as the background value.
  Ignored if there is no \code{status} vector.}
  \item{pch}{vector or list of plotting characters.
  Ignored is there is no \code{status} vector.}
  \item{col}{numeric or character vector of colors, of the same length as \code{values}.
  Defaults to \code{1+1:length(values)}.
  Ignored if there is no \code{status} vector.}
  \item{cex}{numeric vector of plot symbol expansions.
  If a vector, then of the same length as \code{values}. 
  Ignored if there is no \code{status} vector.}
  \item{legend}{character string giving position to place legend.
  See \code{\link{legend}} for possible values.
  Can also be logical, with \code{FALSE} meaning no legend.
  Ignored if there is no \code{status} vector.}
  \item{zero.weights}{logical, should spots with zero or negative weights be plotted?}
  \item{pch.bg}{plotting character for background (non-highlighted) points.}
  \item{col.bg}{color for background (non-highlighted) points.}
  \item{cex.bg}{plot symbol expansion for background (non-highlighted) points.}
  \item{\dots}{The \code{plot} methods pass other arguments to \code{plotWithHighlights}, and \code{plotWithHighlights} passes other arguments to \code{plot.default}.}
}

\details{
An MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values).
If \code{x} is an \code{RGList} or \code{MAList} then this function produces an ordinary within-array MA-plot.
If \code{x} is an \code{MArrayLM} object, then the plot is an fitted model MA-plot in which the estimated coefficient is on the y-axis and the average A-value is on the x-axis.

If \code{x} is a \code{EList} object, then this function produces a between-array MA-plot.
An articifial array is produced by averaging all the arrays other than the array specified.
A mean-difference plot is then producing from the specified array and the artificial array.
Note that this procedure reduces to an ordinary mean-difference plot when there are just two arrays total.

The \code{status} vector is intended to specify the control status of each spot, for example \code{"gene"}, \code{"ratio control"}, \code{"house keeping gene"}, \code{"buffer"} and so on.
The vector is often computed using the function \code{\link{controlStatus}} and a spot-types file.
However the function may be used to highlight any subset of spots.

The \code{status} can be included as the component \code{x$genes$Status} instead of being passed as an argument to \code{plot}.
The arguments \code{values}, \code{pch}, \code{col} and \code{cex} can be included as attributes to \code{status} instead of being passed as arguments to \code{plotMA}.

See \code{\link[graphics]{points}} for possible values for \code{pch}, \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} Volume 43.
\url{http://www.statsci.org/smyth/pubs/limmaPreprint.pdf}
}

\author{Gordon Smyth}

\examples{
A <- runif(1000,4,16)
y <- A + matrix(rnorm(1000*3,sd=0.2),1000,3)
status <- rep(c(0,-1,1),c(950,40,10))
y[,1] <- y[,1] + status
E <- new("EList",list(E=y))
plot(E,array=1,status=status,values=c(-1,1),col=c("blue","red"))
}

\seealso{
\code{\link[limma:plotma]{plotMA}}, \code{\link{plotFB}}, \code{\link{plotMDS}}, \code{\link{plotSA}}

An overview of diagnostic plots available in LIMMA is given in \link{09.Diagnostics}.
}
