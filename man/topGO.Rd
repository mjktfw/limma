\name{topGO}
\alias{topGO}
\title{Table of Top GO Terms}

\description{
Extract top GO terms from goana results.
}

\usage{
topGO(results, ontology = c("BP", "CC", "MF"), sort = NULL, number = 20L)
}

\arguments{
  \item{results}{data frame produced by \code{\link{goana}}.} 
  \item{ontology}{character vector of ontologies to be included in output.  Elements should be one or more of \code{"BP"}, \code{"CC"} or \code{"MF"}.}
  \item{sort}{name of gene set for which results are required.  Should be one of the column names of \code{results}.  Defaults to first set.}
  \item{number}{maximum number of top GO terms to list. For all terms, set \code{number=Inf}.}
}

\details{
This function is organize the output from \code{\link{goana}} into top-tables of the most significant GO terms.
}

\value{
Same as \code{results} but with rows subsetted by Ontology and sorted by the specified p-value.
}

\seealso{
\code{\link{goana}}, \code{\link{goana.MArrayLM}}
}

\author{Gordon Smyth and Yifang Hu}

\examples{
# See goana and goana.MArrayLM examples
}

\keyword{gene set test}
