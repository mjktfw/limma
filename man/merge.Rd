\name{merge.RGList}
\alias{merge.RGList}
\title{Merge RGList Data Objects}
\description{
Merge two microarray data sets represented by RGLists.
}
\usage{
\method{merge}{RGList}(x,y,\dots)
}
\arguments{
  \item{x}{\code{\link{RGList-class}} object with list components \code{R}, \code{G}, \code{Rb} and \code{Gb} containing the foreground and background intensities for each spot on each array.}
  \item{y}{\code{RGList} object, corresponding to the same genes as for \code{x}, possibly in a different order, but with different arrays.}
  \item{\dots}{other arguments are accepted but not used at present}
}
\details{
An \code{RGList} is a list object containing numeric matrices all of the same dimensions.
The RGLists are merged by merging each of the components by row names.
Unlike when using \code{\link{cbind}}, row names are not required to be in the same order or to be unique.
In the case of repeated row names, the order of the rows with repeated names in preserved.
This means that the first occurrence of each name in \code{x$R} is matched with the first occurrence of the same name in \code{y$R}, the second with the second, and so on.
The final vector of row names is the same as in \code{x}.
}
\value{
An \code{RGList} with the same components as \code{x}.
Component matrices have the same rows names as in \code{x} but columns from \code{y} as well as \code{x}.
}
\author{Gordon Smyth}
\seealso{
R base provides a \code{\link[base]{merge}} method for merging data.frames.

An overview of limma commands for reading, subsetting and merging data is given in \link{3.ReadingData}.
}
\examples{
R <- G <- matrix(11:14,4,2)
rownames(R) <- rownames(G) <- c("a","a","b","c")
RG1 <- new("RGList",list(R=R,G=G))

R <- G <- matrix(21:24,4,2)
rownames(R) <- rownames(G) <- c("b","a","a","c")
RG2 <- new("RGList",list(R=R,G=G))

merge(RG1,RG2)
merge(RG2,RG1)
}
\keyword{manip}
