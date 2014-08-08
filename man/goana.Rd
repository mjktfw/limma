\name{goana}
\alias{goana}
\alias{goana.default}
\alias{goana.MArrayLM}
\title{Gene Ontology Analysis of Differentially Expressed Genes}
\description{
Fisher's exact test or Wallenius' noncentral hypergeometric test on differentially expressed genes for over-representation of gene ontology (GO) terms.
}
\usage{
\method{goana}{default}(de, universe = NULL, species = "Hs", weights = NULL, \dots)
\method{goana}{MArrayLM}(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, 
      species = "Hs", trend = FALSE, \dots)
}
\arguments{
  \item{de}{Differentially expressed genes.  Can be a \code{MArrayLM} fit object, or a vector or list of Entrez Gene IDs of differentially expressed genes.}
  \item{universe}{vector specifying the Entrez Gene identifiers of the universe. If \code{NULL}, all Entrez Gene IDs with any gene ontology annotation will be the universe.}
  \item{species}{character string specifying the species. Possible values are \code{"Hs"}, \code{"Mm"}, \code{"Rn"} or \code{"Dm"}.}
  \item{weights}{probability weighting of the universe.}
  \item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest.}
  \item{geneid}{Entrez Gene identifiers. Either a vector of length \code{nrow(de)} or the name of the column of \code{de$genes} containing the Entrez Gene IDs.}
  \item{FDR}{numeric. False discovery rate of differentially expressed genes less than this cutoff.}
  \item{trend}{adjust analysis for gene length or abundance?
  Can be logical, or a numeric vector of covariate values, or the name of the column of \code{de$genes} containing the covariate values.
  If \code{TRUE}, then \code{de$Amean} is used as the covariate.}
  \item{\dots}{further arguments.}
}
\details{
Performs a Gene Ontology enrichment analysis for one for more gene lists using the appropriate Bioconductor organism package.
The gene lists must be supplied as Entrez Gene IDs.

If applied to a linear model fit, then the differentially expressed genes are extracted automatically from the fit object, and analyses are done separately for the up and down significant genes.

If \code{trend=FALSE}, the function computes one-sided hypergeometric tests.
If \code{trend=TRUE} or a covariate is supplied, then a trend is fitted to the differential expression results and the method of Young et al (2010) is used to adjust for this trend.
The adjusted test uses Wallenius' noncentral hypergeometric distribution.
}
\value{
  A dataframe with rows for GO IDs and the following possible columns:
  \item{Term}{GO terms.}
  \item{Ont}{GO term ontology from BP, CC and MF.}
  \item{N}{Number of genes in the GO term.}
  \item{Up}{Number of up regulated differentially expressed genes for testing.}
  \item{Down}{Number of down regulated differentially expressed genes for testing.}
  \item{DE1}{Number of differentially expressed genes for testing.}
  \item{P.Up}{P value from hypergeometric test or Wallenius' noncentral hypergeometric test of up regulated differentially expressed genes for over representation of GO terms.}
  \item{P.Down}{P value from hypergeometric test or Wallenius' noncentral hypergeometric test of down regulated differentially expressed genes for over representation of GO terms.}
  \item{P.DE1}{P value from hypergeometric test or Wallenius' noncentral hypergeometric test of differentially expressed genes for over representation of GO terms.}
}

\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{topGO}}

The goseq package implements a similar GO analysis.
The goseq version will work with a variety of gene identifiers, not only Entrez Gene as here, and includes a database of gene length information for various species.

The gostats package also does GO analyses with some different options.
}
\author{Gordon Smyth and Yifang Hu}
\examples{

\dontrun{

fit <- lmFit(y, design)
fit <- eBayes(fit)

# goana without adjusting for gene length bias
go.fisher <- goana(fit)
topGO(go.fisher, sort = "up")
topGO(go.fisher, sort = "down")

# goana adjusting for gene length bias 
go.len <- goana(fit, geneid = "GeneID", trend = "Length")
topGO(go.len, sort = "up")
topGO(go.len, sort = "down")

# goana adjusting for gene abundance
go.abund <- goana(fit, geneid = "GeneID", trend = TRUE)
topGO(go.abund, sort = "up")
topGO(go.abund, sort = "down")

# goana.default
go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.de, ontology = "BP", sort = "DE3")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")

}

}

\keyword{gene set test}
