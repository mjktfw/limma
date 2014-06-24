\name{goana}
\alias{goana}
\title{Gene Ontology Analysis of Differentially Expressed Genes}
\description{
Hypergeometric tests on up and down differentially expressed genes for over-representation of GO Terms.
}
\usage{
goana(fit, coef = ncol(fit), geneid = "GeneID", FDR = 0.05, species = "Hs")
}
\arguments{
  \item{fit}{should be an object of class \code{MArrayLM} as produced by \code{lmFit} and \code{eBayes}.}
  \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest.}
  \item{geneid}{Entrez gene identifiers. Either a vector of length nrow(fit) or the name of the column of fit$genes containing the Entrez Gene IDs.}
  \item{FDR}{numeric. False discovery rate of differentially expressed genes less than this cutoff.}
  \item{species}{character string specifying the species. Possible values are \code{"Hs"}, \code{"Mm"}, \code{"Rn"} or \code{"Dm"}.}
}
\details{
This function is used in conjunction with \code{\link{lmFit}}, \code{\link{eBayes}} and \code{\link{topGO}}.
}
\value{
  A dataframe with rows for GO IDs and the following columns:
  \item{Term}{GO terms.}
  \item{Ont}{GO term ontology from BP, CC and MF.}
  \item{N}{Number of genes in the GO term.}
  \item{Up}{Number of up regulated differentially expressed genes for testing.}
  \item{Down}{Number of down regulated differentially expressed genes for testing.}
  \item{P.Up}{P value from hypergeometric test of up regulated differentially expressed genes for over representation of GO terms.}
  \item{P.Down}{P value from hypergeometric test of down regulated differentially expressed genes for over representation of GO terms.}
}

\seealso{
\code{\link{topGO}}

The gostats package also does GO analyses with some extra options.
}
\author{Gordon Smyth and Yifang Hu}
\examples{

\dontrun{

fit <- lmFit(y,design)
fit <- eBayes(fit)
go.results <- goana(fit)
topGO(go.results, sort = "up")
topGO(go.results, sort = "down")

}

}

\keyword{gene set test}
