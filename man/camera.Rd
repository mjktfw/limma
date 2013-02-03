\name{camera}
\alias{camera}
\alias{camera.EList}
\alias{camera.MAList}
\alias{camera.default}
\alias{interGeneCorrelation}
\title{Competitive Gene Set Test Accounting for Inter-gene Correlation}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, accounting for inter-gene correlation.
}
\usage{
\method{camera}{EList}(y, index, design, contrast=ncol(design), weights=NULL, use.ranks=FALSE, allow.neg.cor=TRUE, trend.var=FALSE, sort=TRUE)
interGeneCorrelation(y, design)
}
\arguments{
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any data object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
  Rows correspond to probes and columns to samples.}
  \item{index}{an index vector or a list of index vectors.  Can be any vector such that \code{y[index,]} selects the rows corresponding to the test set.}
  \item{design}{design matrix.}
  \item{contrast}{contrast of the linear model coefficients for which the test is required. Can be an integer specifying a column of \code{design}, or else a numeric vector of same length as the number of columns of \code{design}.}
  \item{weights}{can be a numeric matrix of individual weights, of same size as \code{y}, or a numeric vector of array weights with length equal to \code{ncol(y)}.}
  \item{use.ranks}{do a rank-based test (\code{TRUE}) or a parametric test (\code{FALSE}?}
  \item{allow.neg.cor}{should reduced variance inflation factors be allowed for negative correlations?}
  \item{trend.var}{logical, should an empirical Bayes trend be estimated?  See \code{\link{eBayes}} for details.}
  \item{sort}{logical, should the results be sorted by p-value?}
}

\details{
\code{camera} and \code{interGeneCorrelation} implement methods proposed by Wu and Smyth (2012).
\code{camera} performs a \emph{competitive} test in the sense defined by Goeman and Buhlmann (2007).
It tests whether the genes in the set are highly ranked in terms of differential expression relative to genes not in the set.
It has similar aims to \code{geneSetTest} but accounts for inter-gene correlation.
See \code{\link{roast}} for an analogous \emph{self-contained} gene set test.

The function can be used for any microarray experiment which can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model by giving the vector of test statistic values.

\code{camera} estimates p-values after adjusting the variance of test statistics by an estimated variance inflation factor.
The inflation factor depends on estimated genewise correlation and the number of genes in the gene set.

\code{interGeneCorrelation} estimates the mean pair-wise correlation between a set of genes.
}

\value{
\code{camera} returns a data.frame with a row for each set and the following columns:
\item{NGenes}{number of genes in set}
\item{Correlation}{inter-gene correlation}
\item{Direction}{direction of change (\code{"Up"} or \code{"Down"})}
\item{PValue}{two-tailed p-value}
\item{FDR}{Benjamini and Hochberg FDR adjusted p-value}

\code{interGeneCorrelation} returns a list with components:
\item{vif}{variance inflation factor}
\item{correlation}{inter-gene correlation}
}

\author{Di Wu and Gordon Smyth}

\references{
Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\url{http://nar.oxfordjournals.org/content/40/17/e133}

Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 
}

\seealso{
\code{\link{rankSumTestWithCorrelation}},
\code{\link{geneSetTest}},
\code{\link{roast}},
\code{\link{romer}}.
}

\examples{
y <- matrix(rnorm(1000*6),1000,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

# First set of 20 genes are genuinely differentially expressed
index1 <- 1:20
y[index1,4:6] <- y[index1,4:6]+1

# Second set of 20 genes are not DE
index2 <- 21:40
 
camera(y, index1, design)
camera(y, index2, design)

camera(y, list(set1=index1,set2=index2), design)
}
\keyword{htest}

