\name{QualityWeights}
\alias{QualityWeights}
\alias{wtarea}
\alias{wtflags}
\alias{wtIgnore.Filter}
\title{Spot Quality Weights}
\description{
Functions to calculate quality weights for individual spots based on image analyis output file.
}
\usage{
wtarea(ideal=c(160,170))
wtflags(flagged.wt=0.1)
wtIgnore.Filter
}
\arguments{
  \item{ideal}{numeric vector giving the ideal value or range of values}
  \item{flagged.wt}{weight to be given to flagged spots}
}
\details{
These functions can be passed as an argument to \code{read.maimages} to construct quality weights as the microarray data is read in.

\code{wtarea} downweights unusually small or large spots and is designed for SPOT output.
It gives weight 1 to spots which have areas in the ideal range, given in pixels, and linearly downweights spots which are smaller or larger than this range.

\code{wtflags} is designed for GenePix output and gives weight \code{flagged.wt} to spots with \code{Flags} value less than zero. 

\code{wtIgnore.Filter} is designed for QuantArray output and gives weight equal to the column \code{Ignore Filter} produced by QuantArray.
}
\value{
A function which takes a dataframe or matrix as argument and produces a numeric vector of weights between 0 and 1
}
\author{Gordon Smyth}
\seealso{
\code{\link{read.maimages}}, \code{\link{rg.series.spot}}
}
\examples{
#  Read in spot output files from current directory and give full weight to 165 pixel spots
#  RG <- read.maimages(fnames,source="spot",wt.fun=wtarea(165))
#  normalize will downweight spots according to weights found in RG
#  MA <- normalize(RG,layout)
}
\keyword{regression}
