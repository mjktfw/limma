\name{read.maimages}
\alias{read.maimages}
\title{Read RG list from Image Analysis Output Files}
\description{
Extracts an RG list from a series of microarray image analysis output files
}
\usage{
read.maimages(files,source="spot",path=NULL,ext=NULL,names=files,columns=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote="\"",...)

}
\arguments{
  \item{files}{character vector giving the names of the files containing the image analysis output}
  \item{source}{the image analysis program that produced the output files.  Choices are \code{"spot"}, \code{"spot.close.open"}, \code{"genepix"} or \code{"quantarray"}.}
  \item{path}{character string giving the directory containing the files. Can be omitted if the files are in the current working directory.}
  \item{ext}{optional extension to be added to each file name}
  \item{names}{names to be associated with each array as column name}
  \item{columns}{list with fields \code{Rf}, \code{Gf}, \code{Rb} and \code{Gb} giving the column names to be used for red and green foreground and background.  This is not usually specified by the user, if it is, then it over-rides \code{source}.}
  \item{wt.fun}{function to calculate quality weights}
  \item{verbose}{\code{TRUE} to report each time a file in read in}
  \item{sep}{the field separator character}
  \item{quote}{the set of quotting characters}
  \item{\dots}{any other arguments to be passed to \code{read.table}}
}
\details{
This function extracts the foreground and background intensities from a series of files, produced by an image analysis program, and assembles them in the components of one list.

Spot quality weights may also be extracted using a ready-made or a user-supplied weight function.
See \code{\link{QualityWeights}} for provided weight functions.
}
\value{
A list containing the components
  \item{R}{matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{matrix containing the green channel background intensities for each spot for each array.}
  \item{weights}{spot quality weights, if \code{wt.fun} is given}
}
\author{Gordon Smyth}
\seealso{
\code{\link{read.series}}, \code{\link{rg.spot}}, \code{\link{QualityWeights}}, \code{\link[base]{read.table}}
}
\examples{
#  Read all .gpr files from current working directory
#
#  files <- dir(pattern="*\\.gpr$")
#  RG <- read.maimages(files,"genepix")
}
\keyword{file}
