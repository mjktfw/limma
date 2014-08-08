\title{Differential splicing plot}
\name{plotSplice}
\alias{plotSplice}
\description{
Plot relative log-fold changes by exons for the specified gene.
}
\usage{
plotSplice(fit, coef=ncol(fit), geneid=NULL, genecolname=NULL, rank=1L, FDR = 0.05)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is assessed.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{genecolname}{column name of \code{fit$genes} containing gene IDs. Defaults to \code{fit$genecolname}.}
  \item{rank}{integer, if \code{geneid=NULL} then this ranked gene will be plotted.}
  \item{FDR}{numeric, mark exons with false discovery rate less than this cutoff.}
}

\details{
Plot relative log2-fold-changes by exon for the specified gene.
The relative logFC is the difference between the exon's logFC and the overall logFC for the gene, as computed by \code{diffSplice}.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth and Yifang Hu}
\seealso{
\code{\link{diffSplice}}

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\examples{# See diffSplice}
\keyword{rna-seq}
