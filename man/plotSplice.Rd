\title{Plot exons on differentially spliced gene}
\name{plotSplice}
\alias{plotSplice}
\description{
Plot exons of differentially spliced gene.
}
\usage{
plotSplice(fit, coef=ncol(fit), geneid=NULL, rank=1L, FDR = 0.05)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is assessed.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{rank}{integer, if \code{geneid=NULL} then this ranked gene will be plotted.}
  \item{FDR}{numeric, mark exons with false discovery rate less than this cutoff.}
}

\details{
Plots interaction log2-fold-change by exon for the specified gene.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth and Yifang Hu}
\seealso{
\code{\link{diffSplice}}

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\examples{# See diffSplice}
