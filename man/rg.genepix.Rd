\name{rg.genepix}
\alias{rg.genepix}
\title{Extract RG list from Genepix output files}
\description{
Extracts an RG list from a series of Genepix output files
}
\usage{
rg.quantarray(slides,names.slides=names(slides),suffix="gpr")
}
\arguments{
  \item{slides}{Character vector giving the names of the data frames containing the Spot output.}
  \item{names.slides}{Names to be associated with each slide as column name.}
  \item{suffix}{The dataframe names are assumed to have this suffix added to names in \code{slides}.}
}
\details{
This function extracts the foreground and background intensities from a series of data frames and assembles them in the components of one list.
}
\value{
A list containing the components
  \item{R}{A matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{A matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{A matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{A matrix containing the green channel background intensities for each spot for each array.}
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\keyword{IO}
