\name{venn}
\alias{vennCounts}
\alias{vennDiagram}
\title{Venn Diagrams}
\description{
Compute classification counts or plot classification counts in a Venn diagram.
}
\usage{
vennCounts(classification)
vennDiagram(object, names, \dots)
}
\arguments{
  \item{classification}{classification matrix of 0's and 1's indicating significance of a test.
  Usually created by \code{\link{classifyTests}}.}
  \item{object}{either a classification matrix or a \code{VennCounts} object produced by \code{vennCounts}.}
  \item{names}{optional character vector giving names for the sets or contrasts}
  \item{\dots}{any other arguments are passed to \code{plot}}
}
\value{
\code{vennCounts} produces a \code{VennCounts} object, which is a numeric matrix with last column \code{"Counts"} giving counts for each possible vector outcome.
}
\seealso{
An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\author{Gordon Smyth and James Wettenhall}
\examples{
tstat <- matrix(rt(300,df=10),100,3)
tstat[1:33,] <- tstat[1:33,]+2
clas <- classifyTests(tstat,df=10,p.value=0.05)
a <- vennCounts(clas)
print(a)
vennDiagram(a)
}
\keyword{htest}
