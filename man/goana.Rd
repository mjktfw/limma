\name{goana}
\alias{goana}
\alias{goana.default}
\alias{goana.MArrayLM}
\alias{kegga}
\alias{kegga.default}
\alias{kegga.MArrayLM}
\title{Gene Ontology or KEGG Pathway Analysis}

\description{
Test for over-representation of gene ontology (GO) terms or KEGG pathways in one or more sets of genes, optionally adjusting for abundance or gene length bias.
}

\usage{
\method{goana}{MArrayLM}(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, \dots)
\method{goana}{default}(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL,
      plot=FALSE, \dots)
\method{kegga}{MArrayLM}(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, \dots)
\method{kegga}{default}(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL,
      plot=FALSE, \dots)
}

\arguments{
  \item{de}{a vector of Entrez Gene IDs, or a list of such vectors, or an \code{MArrayLM} fit object.}
  \item{coef}{column number or column name specifying for which coefficient or contrast differential expression should be assessed.}
  \item{geneid}{Entrez Gene identifiers. Either a vector of length \code{nrow(de)} or the name of the column of \code{de$genes} containing the Entrez Gene IDs.}
  \item{FDR}{false discovery rate cutoff for differentially expressed genes. Numeric value between 0 and 1.}
  \item{species}{species identifier. Possible values are \code{"Hs"}, \code{"Mm"}, \code{"Rn"} or \code{"Dm"}.}
  \item{trend}{adjust analysis for gene length or abundance?
  Can be logical, or a numeric vector of covariate values, or the name of the column of \code{de$genes} containing the covariate values.
  If \code{TRUE}, then \code{de$Amean} is used as the covariate.}
  \item{universe}{vector specifying the set of Entrez Gene identifiers to be the background universe.
  If \code{NULL} then all Entrez Gene IDs associated with any gene ontology term will be used as the universe.}
  \item{prior.prob}{optional numeric vector of the same length as \code{universe} giving the prior probability that each gene in the universe appears in a gene set. Will be computed from \code{covariate} if the latter is provided.}
  \item{covariate}{optional numeric vector of the same length as \code{universe} giving a covariate against which \code{prior.prob} should be computed.}
  \item{plot}{logical, should the \code{prior.prob} vs \code{covariate} trend be plotted?}
  \item{\dots}{any other arguments in a call to the \code{MArrayLM} method are passed to the default method.}
}

\details{
These functions performs a over-representation analysis for Gene Ontology terms or KEGG pathways in a list of Entrez Gene IDs.
The default method accepts the gene list as a vector of gene IDs,
while the \code{MArrayLM} method extracts the gene lists automatically from a linear model fit object.

\code{goana} uses annotation from the appropriate Bioconductor organism package.
\code{kegga} uses the KEGGREST package to access KEGG pathway annotation.

The default method accepts a vector \code{prior.prob} giving the prior probability that each gene in the universe appears in a gene set.
This vector can be used to correct for unwanted trends in the differential expression analysis associated with gene length, gene abundance or any other covariate.
The \code{MArrayLM} object computes the \code{prior.prob} vector automatically when \code{trend} is non-\code{NULL}.

If \code{prior.prob=NULL}, the function computes one-sided hypergeometric tests equivalent to Fisher's exact test.
If prior probabilities are specified, then a test based on the Wallenius' noncentral hypergeometric distribution is used to adjust for the relative probability that each gene will appear in a gene set, following the approach of Young et al (2010).

The \code{MArrayLM} methods performs over-representation analyses for the up and down differentially expressed genes from a linear model analysis.
In this case, the universe is all the genes found in the fit object.

\code{trend=FALSE} is equivalent to \code{prior.prob=NULL}.
If \code{trend=TRUE} or a covariate is supplied, then a trend is fitted to the differential expression results and this is used to set \code{prior.prob}.

The statistical approach provided here is the same as that provided by the goseq package, with one methodological difference and a few restrictions.
Unlike the goseq package, the gene identifiers here must be Entrez Gene IDs and the user is assumed to be able to supply gene lengths if necessary.
The goseq package has additional functionality to convert gene identifiers and to provide gene lengths.
The only methodological difference is that \code{goana} and \code{kegga} computes gene length or abundance bias using \code{tricubeMovingAverage} instead of monotonic regression.
While \code{tricubeMovingAverage} does not enforce monotonicity, it has the advantage of numerical stability when \code{de} contains only a small number of genes.
}

\value{
The \code{goana} default method produces a data frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{number of genes in the GO term.}
  \item{DE}{number of genes in the \code{DE} set.}
  \item{P.DE}{p-value for over-representation of the GO term in the set.}
The last two column names above assume one gene set with the name \code{DE}.
In general, there will be a pair of such columns for each gene set and the name of the set will appear in place of \code{"DE"}.

The \code{goana} method for \code{MArrayLM} objects produces a data frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{number of genes in the GO term.}
  \item{Up}{number of up-regulated differentially expressed genes.}
  \item{Down}{number of down-regulated differentially expressed genes.}
  \item{P.Up}{p-value for over-representation of GO term in up-regulated genes.}
  \item{P.Down}{p-value for over-representation of GO term in down-regulated genes.}

The row names of the data frame give the GO term IDs.

The output from \code{kegga} is the same except that row names become KEGG pathway IDs, \code{Term} becomes \code{Pathway} and there is no \code{Ont} column.
}

\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{topGO}}, \code{\link{topKEGG}}

The goseq package provides an alternative implementation of methods from Young et al (2010).
Unlike the limma functions documented here, goseq will work with a variety of gene identifiers and includes a database of gene length information for various species.

The gostats package also does GO analyses without adjustment for bias but with some other options.

See \link{10.GeneSetTests} for a description of other functions used for gene set testing.
}

\author{Gordon Smyth and Yifang Hu}

\examples{
\dontrun{
## Linear model usage:

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

# GO analysis adjusting for gene length bias
# (assuming that y$genes$Length contains gene lengths)
go.len <- goana(fit, geneid = "GeneID", trend = "Length")
topGO(go.len, sort = "up")
topGO(go.len, sort = "down")

## Default usage with a gene list:

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.de, ontology = "BP", sort = "DE3")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")
}
}

\keyword{gene set test}
