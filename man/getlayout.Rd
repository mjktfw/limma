\name{getLayout}
\alias{getLayout}
\title{Extract the Print Layout of an Array from the GAL File}
\description{
From the Block, Row and Column information in the GAL file, determine the number of tip rows and columns on the array and the number of spot rows and columns within each tip group.
}
\usage{
getLayout(gal)
}
\arguments{
  \item{gal}{data.frame containing the GAL, i.e., giving the position and gene identifier of each spot}
}
\details{
A GAL file is a list of genes and associated information produced by an Axon microarray scanner.
This function assumes that the data.frame contains columns \code{Block}, \code{Column} and \code{Row}.
The number of tip columns is not determinable from the GAL but is assumed to be four.
}
\value{
A list with components
  \item{ngrid.r}{number of tip rows on the print head}
  \item{ngrid.c}{number of tip columns on the print head}
  \item{nspot.r}{number of rows of spots printed with each tip}
  \item{nspot.c}{number of columns of spots printed with each tip}
}
\author{Gordon Smyth}
\seealso{
\code{\link[marrayTools]{gpTools}}.

An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\examples{
# gal <- readGAL()
# layout <- getLayout(gal)
}
\keyword{IO}
