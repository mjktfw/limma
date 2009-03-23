\name{cbind}
\alias{cbind.RGList}
\alias{cbind.MAList}
\alias{cbind.EList}
\alias{rbind.RGList}
\alias{rbind.MAList}
\alias{rbind.EList}
\title{Combine RGList, MAList or EList Objects}
\description{
Combine a set of \code{RGList}, \code{MAList} or \code{EList} objects.
}
\usage{
\method{cbind}{RGList}(\dots, deparse.level=1)
\method{rbind}{RGList}(\dots, deparse.level=1)
}
\arguments{
  \item{\dots}{\code{RGList}, \code{MAList} or \code{EList} objects.}
  \item{deparse.level}{not currently used, see \code{\link[base]{cbind}} in the base package}
}
\details{
\code{cbind} combines data objects assuming the same probes in the same order but different arrays.
\code{rbind} combines data objects assuming equivalent arrays, i.e., the same RNA targets, but different probes.

For \code{cbind}, the matrices of expression data from the individual objects are cbinded.
The data.frames of target information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
For \code{rbind}, the matrices of expression data are rbinded while the target information, in any, is unchanged.
}
\value{
An \code{\link[limma:rglist]{RGList}}, \code{\link[limma:malist]{MAList}} or \code{\link[limma:EList]{EList}} object holding data from all the arrays and all genes from the individual objects.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{cbind}} in the base package.
  
  \link{03.ReadingData} gives an overview of data input and manipulation functions in LIMMA.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A1","A2")
MA1 <- new("MAList",list(M=M,A=A))

M <- A <- matrix(21:24,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("B1","B2")
MA2 <- new("MAList",list(M=M,A=A))

cbind(MA1,MA2)
}
\keyword{manip}
