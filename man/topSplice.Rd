\title{Top table of differentially spliced genes or exons}
\name{topSplice}
\alias{topSplice}
\description{
Top table ranking the most differentially spliced genes or exons.
}
\usage{
topSplice(fit, coef=ncol(fit), level="exon", number=10, FDR=1)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is 
assessed.}
  \item{level}{character string, should the table be by \code{"exon"} or by 
\code{"gene"}.}
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this 
cutoff.}
}

\details{
Ranks genes by the Plots interaction log2-fold-change by exon for the specified gene.
}

\value{A data.frame with any annotation columns found in \code{fit} plus the following 
columns
  \item{logFC}{log2-fold change of exon vs other exons for the same gene (if 
\code{level="exon"})}
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
