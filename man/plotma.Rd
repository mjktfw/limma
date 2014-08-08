\title{MA-Plot}
\name{plotMA}
\alias{plotMA}
\alias{plotMA.RGList}
\alias{plotMA.MAList}
\alias{plotMA.EList}
\alias{plotMA.MArrayLM}
\alias{plotMA.default}
\description{
Creates an MA-plot with color coding for control spots.
}
\usage{
\method{plotMA}{default}(MA, array = 1, xlab = "A", ylab = "M", main = colnames(MA)[array],
       xlim = NULL, ylim = NULL, status, values, pch, col, cex, legend = TRUE,
       zero.weights = FALSE, ...)
\method{plotMA}{MArrayLM}(MA, coef = ncol(MA), xlab = "AveExpr", ylab = "logFC",
       main = colnames(MA)[coef], xlim = NULL, ylim = NULL, status, values, pch, col,
       cex, legend = TRUE, zero.weights = FALSE, ...)
}
\arguments{
  \item{MA}{an \code{RGList}, \code{MAList}, \code{EList}, \code{ExpressionSet} or \code{MArrayLM} object.
  Alternatively a numeric \code{matrix}.}
  \item{array}{integer giving the array to be plotted.}
  \item{coef}{integer giving the linear model coefficient to be plotted.}
  \item{xlab}{character string giving label for x-axis}
  \item{ylab}{character string giving label for y-axis}
  \item{main}{character string giving title for plot}
  \item{xlim}{numeric vector of length 2 giving limits for x-axis, defaults to min and max of the data}
  \item{ylim}{numeric vector of length 2 giving limits for y-axis, defaults to min and max of the data}
  \item{status}{character vector giving the control status of each spot on the array, of same length as the number of rows of \code{MA$M}.
  If omitted, all points are plotted in the default color, symbol and size.}
  \item{values}{character vector giving values of \code{status} to be highlighted on the plot. Defaults to unique values of \code{status}.
  Ignored if there is no \code{status} vector.}
  \item{pch}{vector or list of plotting characters. Default is integer code 16 which gives a solid circle.
  Ignored is there is no \code{status} vector.}
  \item{col}{numeric or character vector of colors, of the same length as \code{values}. Defaults to \code{1:length(values)}.
  Ignored if there is no \code{status} vector.}
  \item{cex}{numeric vector of plot symbol expansions, of the the same length as \code{values}. 
  Defaults to 0.3 for the most common status value and 1 for the others.
  Ignored if there is no \code{status} vector.}
  \item{legend}{logical, should a legend of plotting symbols and colors be included. Can also be a character string giving position to place legend. Ignored if there is no \code{status} vector.}
  \item{zero.weights}{logical, should spots with zero or negative weights be plotted?}
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
An MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values).
For two color data objects, a within-array MA-plot is produced with the M and A values computed from the two channels for the specified array.
This is the same as a mean-difference plot (\code{\link{mdplot}}) with the red and green log2-intensities of the array providing the two columns.

For single channel data objects, then a between-array MA-plot is produced.
An articifial array is produced by averaging all the arrays other than the array specified.
A mean-difference plot is then producing from the specified array and the artificial array.
Note that this procedure reduces to an ordinary mean-difference plot when there are just two arrays total.

If \code{MA} is an \code{MArrayLM} object, then the plot is an fitted model MA-plot in which the estimated coefficient is on the y-axis and the average A-value is on the x-axis.


The \code{status} vector is intended to specify the control status of each spot, for example "gene", "ratio control", "house keeping gene", "buffer" and so on.
The vector is often computed using the function \code{\link{controlStatus}} and a spot-types file.
However the function may be used to highlight any subset of spots.

The \code{status} can be included as the component \code{MA$genes$Status} instead of being passed as an argument to \code{plotMA}.
The arguments \code{values}, \code{pch}, \code{col} and \code{cex} can be included as attributes to \code{status} instead of being passed as arguments to \code{plotMA}.

See \code{\link[graphics]{points}} for possible values for \code{pch}, \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}
\references{See \url{http://www.statsci.org/micrarra/refs/maplots.html}}
\author{Gordon Smyth}
\examples{
MA <- new("MAList")
MA$A <- runif(300,4,16)
MA$M <- rt(300,df=3)
status <- rep("Gene",300)
status[1:3] <- "M=0"
MA$M[1:3] <- 0
status[4:6] <- "M=3"
MA$M[4:6] <- 3
status[7:9] <- "M=-3"
MA$M[7:9] <- -3
plotMA(MA,main="MA-Plot with Simulated Data",
       status=status,values=c("M=0","M=3","M=-3"),col=c("blue","red","green"))

#  Same as above
attr(status,"values") <- c("M=0","M=3","M=-3")
attr(status,"col") <- c("blue","red","green")
plotMA(MA,main="MA-Plot with Simulated Data",status=status)

#  Same as above
MA$genes$Status <- status
plotMA(MA,main="MA-Plot with Simulated Data")
}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}
