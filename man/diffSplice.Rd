\name{diffSplice}
\alias{diffSplice}
\title{Test for Differential Splicing}
\description{Given a linear model fit at the exon level, test for differences in exon retention between experimental conditions.}
\usage{
diffSplice(fit, geneid, exonid=NULL, verbose=TRUE)
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit}. Rows should correspond to exons.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{exon identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the exon identifiers.}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
}
\value{
An object of class \code{MArrayLM} containing both exon level and gene level tests.
Results are sorted by geneid and by exonid within gene.
  \item{coefficients}{numeric matrix of coefficients of same dimensions as \code{fit}. Each coefficient is the difference between the log-fold-change for that exon versus the average log-fold-change for all other exons for the same gene.}
  \item{t}{numeric matrix of moderated t-statistics}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{genes}{data.frame of exon annotation}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs}
  \item{gene.F}{numeric vector of moderated F-statistics, one for each gene.}
  \item{gene.F.p.value}{numeric vector of p-values corresponding to \code{gene.F}}
}
\details{
This function is often used in conjunction with voom.
}
\seealso{
\code{\link{voom}}
}
\author{Gordon Smyth and Charity Law}

\examples{
\dontrun{
v <- voom(dge,design)
fit <- lmFit(v,design)
ex <- diffSplice(fit,geneid="EntrezID")
topSplice(ex)
plotSplice(ex)
}
}
