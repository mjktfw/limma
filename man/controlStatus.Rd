\name{controlStatus}
\alias{controlStatus}
\title{Get Spot Status from Spot Types}
\description{
Determine the type (or status) of each spot in the gene list.
}
\usage{
controlStatus(types, genes)
}
\arguments{
  \item{types}{dataframe containing spot type specifiers, usually input using \code{readSpotTypes}}
  \item{genes}{dataframe containing gene IDs and Names, or an \code{RGList} or \code{MAList} containing this dataframe}
}
\details{
This function matches up the regular expressions associated with spot types with the gene list.
}
\value{
Character vector specifying the type (or status) of each spot on the array
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\keyword{IO}
