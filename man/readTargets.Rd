\name{readTargets}
\alias{readTargets}
\title{Read RNA Targets File}
\description{
Read a RNA targets file into dataframe.
}
\usage{
readTargets(file="Targets.txt", path=NULL, sep="\t")
}
\arguments{
  \item{file}{character string giving the name of the targets file.}
  \item{path}{character string giving the directory containing the file.
       Can be omitted if the file is in the current working irectory.} 
  \item{sep}{the field separator character}
}
\details{
The targets file is a text file with rows corresponding to microarrays and columns \code{Cy3} and \code{Cy5} specifying which RNA samples are hybridized to which channel of each microarray.
Other columns are optional.
}
\value{
A data frame including columns
  \item{Cy3}{character vector giving names of RNA samples}
  \item{Cy5}{character vector giving names of RNA samples}
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\keyword{IO}
