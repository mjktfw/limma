\title{Image Plot of Microarray Statistics}
\name{imageplot}
\alias{imageplot}
\description{
Creates an image of shades
of gray or colors, that represents the values of a statistic for each
spot on the array. The statistic can be a log intensity ratio, quality
information such as spot size or shape, or a t-statistic. This function
can be used to explore whether there are any spatial effects in the data.}
\usage{
imageplot(z, layout = list(ngrid.r = 12, ngrid.c = 4, nspot.r = 26, 
nspot.c = 26), low = NULL, high = NULL, ncolors = 123, zerocenter = NULL, 
zlim = NULL, ...)
}
\arguments{
  \item{z}{a numeric vector or array. This vector can contain any spot 
statistics, such
as log intensity ratios, spot sizes or shapes, or t-statistics. Missing values 
are allowed and will result in blank spots on the image.}
  \item{layout}{a list specifying the dimensions of the spot matrix
and the grid matrix.}
  \item{low}{color associated with low values of \code{z}. May be specified as a character string 
such as \code{"green"}, \code{"white"} etc, or as a rgb vector in which \code{c(1,0,0)} is red, 
\code{c(0,1,0)} is green and \code{c(0,0,1)} is blue. The default value is \code{"green"} if \code{zerocenter=T} or \code{"white"} if \code{zerocenter=F}.}
  \item{high}{color associated with high values of \code{z}. The default value is \code{"red"} if \code{zerocenter=T} or \code{"blue"} if \code{zerocenter=F}.}
  \item{ncolors}{number of color shades used in the image including low and high.}
  \item{zerocenter}{should zero values of \code{z} correspond to a shade exactly halfway between the colors 
low and high? The default is TRUE if \code{z} takes positive and negative values, 
otherwise FALSE.}
  \item{zlim}{numerical vector of length 2 giving the extreme values of \code{z} to associate with 
colors \code{low} and \code{high}. By default \code{zlim} is the range of \code{z}. Any values of \code{z} outside 
the interval \code{zlim} will be truncated to the relevant limit.}
\item{...}{any other arguments will be passed to the function image}
}
\details{
The image follows the layout of an actual microarray slide with the bottom left corner representing the spot (1,1,1,1). 
This function is very similar to the \code{sma} function \code{plot.spatial} but is intended 
to display spatial patterns and artefacts rather than to highlight extreme 
values. The function differs from \code{plot.spatial} most noticeably in that all the 
spots are plotted and the image is plotted from bottom left rather than from 
top left.
}
\value{An image is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{\code{\link[marrayPlots]{maImage}}, \code{\link[base]{image}}.}
\examples{
M <- rnorm(8*4*16*16)
imageplot(M,layout=list(ngrid.r=8,ngrid.c=4,nspot.r=16,nspot.c=16))
}
\keyword{hplot}
