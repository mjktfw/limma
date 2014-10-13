\name{goana}
\alias{goana}
\alias{goana.default}
\title{Gene Ontology Analysis}

\description{
Test for over-representation of gene ontology (GO) terms in one or more sets of genes.
}

\usage{
\method{goana}{default}(de, universe = NULL, species = "Hs", prior.prob = NULL, \dots)
}

\arguments{
  \item{de}{a vector of Entrez Gene IDs, or a list of such vectors.}
  \item{universe}{vector specifying the set of Entrez Gene identifiers to be the background universe.
  If \code{NULL} then all Entrez Gene IDs associated with any gene ontology term will be used as the universe.}
  \item{species}{species identifier. Possible values are \code{"Hs"}, \code{"Mm"}, \code{"Rn"} or \code{"Dm"}.}
  \item{prior.prob}{numeric vector giving the prior probability that each gene in the universe appears in a gene set.}
  \item{\dots}{other arguments are not currently used.}
}

\details{
\code{goana} is an S3 generic function.
The default method performs a Gene Ontology enrichment analysis for one for more gene lists using the appropriate Bioconductor organism package.
The gene lists must be supplied as Entrez Gene IDs.

If \code{prior.prob=NULL}, the function computes one-sided hypergeometric tests equivalent to Fisher's exact test.

The \code{prior.prob} vector can be used to specify the prior probability that each gene in the universe appears in a gene set.
If prior probabilities are specified, then a test based on the Wallenius' noncentral hypergeometric distribution is used to adjust for the relative probability that each gene will appear in a gene set, following the approach of Young et al (2010).
}

\value{
A data frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{number of genes in the GO term.}
  \item{DE1}{number of genes in the \code{DE1} set.}
  \item{P.DE1}{p-value for over-representation of the GO term in the set.}
The last two column names above assume one gene set with the name \code{DE1}.
In general, there will be a pair of such columns for each gene set and the name of the set will appear in place of \code{"DE1"}.

The row names of the data frame give the GO term IDs.
}

\references{
  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{goana.MArrayLM}}, \code{\link{topGO}}

The goseq package implements a similar GO analysis.
The goseq version will work with a variety of gene identifiers, not only Entrez Gene as here, and includes a database of gene length information for various species.

The gostats package also does GO analyses with some different options.
}

\author{Gordon Smyth and Yifang Hu}

\examples{
\dontrun{

go.de <- goana(list(DE1 = EG.DE1, DE2 = EG.DE2, DE3 = EG.DE3))
topGO(go.de, sort = "DE1")
topGO(go.de, sort = "DE2")
topGO(go.de, ontology = "BP", sort = "DE3")
topGO(go.de, ontology = "CC", sort = "DE3")
topGO(go.de, ontology = "MF", sort = "DE3")
}
}

\keyword{gene set test}
