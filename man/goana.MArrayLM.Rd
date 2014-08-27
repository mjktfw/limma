\name{goana.MarrayLM}
\alias{goana.MArrayLM}
\title{Gene Ontology Analysis of Differentially Expressed Genes}

\description{
Test for over-representation of gene ontology (GO) terms in the up and down differentially expressed genes from a linear model fit.
}

\usage{
\method{goana}{MArrayLM}(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, 
      species = "Hs", trend = FALSE, \dots)
}

\arguments{
  \item{de}{an \code{MArrayLM} fit object.}
  \item{coef}{column number or column name specifying for which coefficient or contrast differential expression should be assessed.}
  \item{geneid}{Entrez Gene identifiers. Either a vector of length \code{nrow(de)} or the name of the column of \code{de$genes} containing the Entrez Gene IDs.}
  \item{FDR}{false discovery rate cutoff for differentially expressed genes. Numeric value between 0 and 1.}
  \item{species}{species identifier. Possible values are \code{"Hs"}, \code{"Mm"}, \code{"Rn"} or \code{"Dm"}.}
  \item{trend}{adjust analysis for gene length or abundance?
  Can be logical, or a numeric vector of covariate values, or the name of the column of \code{de$genes} containing the covariate values.
  If \code{TRUE}, then \code{de$Amean} is used as the covariate.}
  \item{\dots}{any other arguments are passed to \code{goana.default}.}
}

\details{
Performs Gene Ontology enrichment analyses for the up and down differentially expressed genes from a linear model analysis.
The Entrez Gene ID must be supplied for each gene.

If \code{trend=FALSE}, the function computes one-sided hypergeometric tests equivalent to Fisher's exact test.

If \code{trend=TRUE} or a covariate is supplied, then a trend is fitted to the differential expression results and the method of Young et al (2010) is used to adjust for this trend.
The adjusted test uses Wallenius' noncentral hypergeometric distribution.
}

\value{
A data frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{number of genes in the GO term.}
  \item{Up}{number of up-regulated differentially expressed genes.}
  \item{Down}{number of down-regulated differentially expressed genes.}
  \item{P.Up}{p-value for over-representation of GO term in up-regulated genes.}
  \item{P.Down}{p-value for over-representation of GO term in down-regulated genes.}
The row names of the data frame give the GO term IDs.
}

\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{goana.default}}, \code{\link{topGO}}

The goseq package implements a similar GO analysis.
The goseq version will work with a variety of gene identifiers, not only Entrez Gene as here, and includes a database of gene length information for various species.

The gostats package also does GO analyses with some different options.
}
\author{Gordon Smyth and Yifang Hu}
\examples{

\dontrun{

fit <- lmFit(y, design)
fit <- eBayes(fit)

# Standard GO analysis
go.fisher <- goana(fit)
topGO(go.fisher, sort = "up")
topGO(go.fisher, sort = "down")

# GO analysis adjusting for gene abundance
go.abund <- goana(fit, geneid = "GeneID", trend = TRUE)
topGO(go.abund, sort = "up")
topGO(go.abund, sort = "down")

# GO analysis adjusting for gene length bias (assuming that y$genes contains a column of gene lengths)
go.len <- goana(fit, geneid = "GeneID", trend = "Length")
topGO(go.len, sort = "up")
topGO(go.len, sort = "down")
}
}

\keyword{gene set test}
