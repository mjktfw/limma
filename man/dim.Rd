\name{dim}
\alias{dim.RGList}
\alias{dim.MAList}
\alias{dim.MArrayLM}
\title{Retrieve the Dimensions of an RGList, MAList or MArrayLM Object}
\description{
Retreive the number of rows (genes) and columns (arrays) for an RGList, MAList or MArrayLM object.
}
\usage{
\method{dim}{RGList}(x)
}
\arguments{
  \item{x}{an object of class \code{RGList}, \code{MAList} or \code{MArrayLM}}
}
\details{
The matrices of expression data from the individual objects are cbinded.
The data.frames of target information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (arrays).
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{dim}} in the base package.
  
  \link{2.Classes} gives an overview of data classes used in LIMMA.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A1","A2")
MA <- new("MAList",list(M=M,A=A))
dim(M)
}
\keyword{array}
