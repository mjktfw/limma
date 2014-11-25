\title{MA-Plot}
\name{plotMA}
\alias{plotMA}
\alias{plotMA.RGList}
\alias{plotMA.MAList}
\alias{plotMA.EListRaw}
\alias{plotMA.EList}
\alias{plotMA.MArrayLM}
\alias{plotMA.default}

\description{
Creates an MA-plot with color coding for control spots.
}

\usage{
\method{plotMA}{default}(object, array = 1, xlab = "Average log-expression",
       ylab = "Expression log-ratio (this sample vs others)",
       main = colnames(object)[array], status=NULL, \dots)
\method{plotMA}{EList}(object, array = 1, xlab = "Average log-expression",
       ylab = "Expression log-ratio (this sample vs others)",
       main = colnames(object)[array], status=object$genes$Status,
       zero.weights = FALSE, \dots)
\method{plotMA}{RGList}(object, array = 1, xlab = "A", ylab = "M",
       main = colnames(object)[array], status=object$genes$Status,
       zero.weights = FALSE, \dots)
\method{plotMA}{MAList}(object, array = 1, xlab = "A", ylab = "M",
       main = colnames(object)[array], status=object$genes$Status,
       zero.weights = FALSE, \dots)
\method{plotMA}{MArrayLM}(object, coef = ncol(object), xlab = "Average log-expression",
       ylab = "log-fold-change", main = colnames(object)[coef],
       status=object$genes$Status, zero.weights = FALSE, \dots)
}

\arguments{
  \item{object}{an \code{RGList}, \code{MAList}, \code{EList}, \code{ExpressionSet} or \code{MArrayLM} object.
  Alternatively a numeric \code{matrix}.}
  \item{array}{integer giving the array to be plotted.}
  \item{coef}{integer giving the linear model coefficient to be plotted.}
  \item{xlab}{character string giving label for x-axis}
  \item{ylab}{character string giving label for y-axis}
  \item{main}{character string giving title for plot}
  \item{status}{vector giving the control status of each spot on the array, of same length as the number of rows of \code{object}.
  If \code{NULL}, then all points are plotted in the default color, symbol and size.}
  \item{zero.weights}{logical, should spots with zero or negative weights be plotted?}
  \item{\dots}{other arguments are passed to \code{\link{plotWithHighlights}}.}
}

\details{
An MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values).
For two color data objects, a within-array MA-plot is produced with the M and A values computed from the two channels for the specified array.
This is the same as a mean-difference plot (\code{\link{mdplot}}) with the red and green log2-intensities of the array providing the two columns.

For single channel data objects, then a between-array MA-plot is produced.
An articifial array is produced by averaging all the arrays other than the array specified.
A mean-difference plot is then producing from the specified array and the artificial array.
Note that this procedure reduces to an ordinary mean-difference plot when there are just two arrays total.

If \code{object} is an \code{MArrayLM} object, then the plot is an fitted model MA-plot in which the estimated coefficient is on the y-axis and the average A-value is on the x-axis.

The \code{status} vector is intended to specify the control status of each spot, for example \code{"gene"}, \code{"ratio control"}, \code{"house keeping gene"}, \code{"buffer"} and so on.
The vector is often computed using the function \code{\link{controlStatus}} and a spot-types file.
However the function may be used to highlight any subset of spots.

The \code{status} can be included as the component \code{object$genes$Status} instead of being passed as an argument to \code{plotMA}.
The arguments \code{values}, \code{pch}, \code{col} and \code{cex} can be included as attributes to \code{status} instead of being passed as arguments to \code{plotMA}.
}

\value{A plot is created on the current graphics device.}

\references{See \url{http://www.statsci.org/micrarra/refs/maplots.html}}

\author{Gordon Smyth}

\examples{
MA <- new("MAList")
MA$A <- runif(300,4,16)
MA$M <- rt(300,df=3)

# Spike-in values
MA$M[1:3] <- 0
MA$M[4:6] <- 3
MA$M[7:9] <- -3

status <- rep("Gene",300)
status[1:3] <- "M=0"
status[4:6] <- "M=3"
status[7:9] <- "M=-3"
values <- c("M=0","M=3","M=-3")
col <- c("blue","red","green")

plotMA(MA,main="MA-Plot with 12 spiked-in points",
       status=status, values=values, col=col)

#  Same as above but setting graphical parameters as attributes
attr(status,"values") <- values
attr(status,"col") <- col
plotMA(MA, main="MA-Plot with 12 spiked-in points", status=status)

#  Same as above but passing status as part of object
MA$genes$Status <- status
plotMA(MA, main="MA-Plot with 12 spiked-in points")

#  Change settings for background points
MA$genes$Status <- status
plotMA(MA, pch.bg=1, cex.bg=0.5)
}

\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\keyword{hplot}
