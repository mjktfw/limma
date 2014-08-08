\title{Top table of differentially spliced genes or exons}
\name{topSplice}
\alias{topSplice}
\description{
Top table ranking the most differentially spliced genes or exons.
}
\usage{
topSplice(fit, coef=ncol(fit), level="hybrid", number=10, FDR=1)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is assessed.}
  \item{level}{character string, possible values are \code{"gene"}, \code{"exon} or \code{"hybrid}.
    \code{"gene"} gives F-tests for each gene.
    \code{"exon"} gives t-tests for each exon.
    \code{"hybrid"} gives gene level results with p-values derived from the exon-level tests.}
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this cutoff.}
}

\details{
Ranks genes or exons by evidence for differential splicing.
The gene-level test is an F-test for any differences in exon usage between experimental conditions.
The exon-level tests are t-tests for differences between each exon and all other exons for the same gene.

The hybrid testing method processes the exon-level p-values to give an overall call of differential splicing for each gene.
It returns the minimum Simes-adjusted p-values for each gene.
It is likely to be more powerful than the genewise F-tests if only a minority of exons for a gene are differentially spliced.
}

\value{A data.frame with any annotation columns found in \code{fit} plus the following columns
  \item{logFC}{log2-fold change of exon vs other exons for the same gene (if \code{level="exon"})}
  \item{t}{moderated t-statistic (if \code{level="exon"})}
  \item{F}{moderated F-statistic (if \code{level="gene"})}
  \item{P.Value}{p-value}
  \item{FDR}{false discovery rate}
}

\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\examples{# See diffSplice}

\keyword{rna-seq}
