\name{read.maimages}
\alias{read.maimages}
\title{Read RGList from Image Analysis Output Files}
\description{
Reads an RGList from a series of microarray image analysis output files
}
\usage{
read.maimages(files,source="spot",path=NULL,ext=NULL,names=files,columns=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote="\"",\dots)
}
\arguments{
  \item{files}{character vector giving the names of the files containing image analysis output or, for Imagene data, a character matrix of names of files}
  \item{source}{character string specifing the image analysis program which produced the output files.  Choices are \code{"spot"}, \code{"spot.close.open"}, \code{"genepix"}, \code{"quantarray"} or \code{"imagene"}.}
  \item{path}{character string giving the directory containing the files. Can be omitted if the files are in the current working directory.}
  \item{ext}{character string giving optional extension to be added to each file name}
  \item{names}{character vector of names to be associated with each array as column name}
  \item{columns}{list with fields \code{Rf}, \code{Gf}, \code{Rb} and \code{Gb} giving the column names to be used for red and green foreground and background.  This is not usually specified by the user but, if it is, it over-rides \code{source}.}
  \item{wt.fun}{function to calculate quality weights}
  \item{verbose}{logical, \code{TRUE} to report each time a file is read in}
  \item{sep}{the field separator character}
  \item{quote}{character string of characters to be treated as quote marks}
  \item{\dots}{any other arguments are passed to \code{read.table}}
}
\details{
This is the main data input function for the LIMMA package.
It extracts the foreground and background intensities from a series of files, produced by an image analysis program, and assembles them into the components of one list.
At present, the image analysis programs SPOT, GenePix and QuantArray are supported.
Other image analysis programs can potentially be accommodated though a user-specified \code{columns} argument.
In the case of SPOT, two possible background estimators are supported:
if \code{source="spot.close.open"} then background intensities are estimated from \code{morph.close.open} rather than \code{morph}.

Spot quality weights may be extracted from the image analysis files using a ready-made or a user-supplied weight function \code{\link[limma:QualityWeights]{wt.fun}}.
\code{wt.fun} may be any user-supplied function which accepts a data.frame argument and returns a vector of non-negative weights.
The columns of the data.frame are as in the image analysis output files.
See \code{\link{QualityWeights}} for provided weight functions.

For Imagene image data the argument \code{files} should be a matrix with two columns.
The first column should contain the names of the files containing green channel (cy3) data and the second column should contain names of files containing red channel (cy5) data.
}
\value{
An \code{\link[limma:rglist]{RGList}} object containing the components
  \item{R}{matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{matrix containing the green channel background intensities for each spot for each array.}
  \item{weights}{spot quality weights, if \code{wt.fun} is given}
}
\author{Gordon Smyth}
\seealso{
\code{read.maimages} is based on \code{\link[base]{read.table}} in the base package.
\code{\link[marrayInput]{read.marrayRaw}} is the corresponding function in the marrayInput package.

An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\examples{
#  Read all .gpr files from current working directory
#  and give weight 0.1 to spots with negative flags

\dontrun{files <- dir(pattern="*\\\\.gpr$")
RG <- read.maimages(files,"genepix",wt.fun=wtflags(0.1))}

#  Read all .spot files from current working director and down-weight
#  spots smaller or larger than 150 pixels

\dontrun{files <- dir(pattern="*\\\\.spot$")
RG <- read.maimages(files,"spot",wt.fun=wtarea(150))}
}
\keyword{file}
