\name{geneSetTest}
\alias{geneSetTest}
\alias{wilcoxGST}
\title{Mean-rank Gene Set Test}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of a given statistic.
Genes are assumed to be independent.
}
\usage{
geneSetTest(index, statistics, alternative="mixed", type="auto", ranks.only=TRUE, nsim=9999)
wilcoxGST(index, statistics, ...)
}

\arguments{
  \item{index}{index vector for the gene set.  This can be a vector of indices, or a logical vector of the same length as \code{statistics} or, in general, any vector such that \code{statistic[index]} gives the statistic values for the gene set to be tested.}
  \item{statistics}{vector, any genewise statistic by which genes can be ranked.}
  \item{alternative}{character string specifying the alternative hypothesis, must be one of \code{"mixed"}, \code{"either"}, \code{"up"} or \code{"down"}. \code{"two.sided"}, \code{"greater"} and \code{"less"} are also permitted as synonyms for \code{"either"}, \code{"up"} and \code{"down"} respectively.}
  \item{type}{character string specifying whether the statistics are signed (t-like, \code{"t"}) or unsigned (F-like, \code{"f"}) or whether the function should make an educated guess (\code{"auto"}).
  If the statistic is unsigned, then it assume that larger statistics are more significant.}
  \item{ranks.only}{logical, if \code{TRUE} only the ranks of the \code{statistics} are used.} 
  \item{nsim}{number of random samples to take in computing the p-value.
  Not used if \code{ranks.only=TRUE}.}
  \item{\ldots}{other arguments are passed to \code{geneSetTest}.}
}

\value{
numeric value giving the estimated p-value.
}

\details{
\code{wilcoxGST} is a synonym for \code{geneSetTest} with \code{ranks.only=TRUE}.
This test procedure was developed by Michaud et al (2008), who called it \emph{mean-rank gene-set enrichment}.

These functions compute a p-value to test the hypothesis that the (indexed) test set of genes tends to be more highly ranked in terms of some test statistic compared to randomly chosen genes.
The statistic might be any statistic of interest, for example a t-statistic or F-statistic for differential expression.

These functions perform \emph{competitive} tests in the sense that genes in the test set are compared to other genes (Goeman and Buhlmann, 2007).
By contrast, a \emph{self-contained} gene set test such as \code{\link{roast}} tests for differential expression for the test set only without regard to other genes on the array.
Like all gene set tests, these functions can be used to detect differential expression for a group of genes, even when the effects are too small or there is too little data to detect the genes individually.
The also provides a means to compare the results between different experiments.

Because it is based on permuting genes, \code{geneSetTest} assumes that the different genes (or probes) are independent.
(Strictly speaking, it assumes that the genes in the set are no more correlated on average than randomly chosen genes.)
This assumption may be reasonable if the gene set is relatively small and if there is relatively little genotypic variation in the data, for example if the data is obtained from genetically identical inbred mice.
The independence assumption may be misleading if the gene set is large or if the data contains a lot of genotypic variation, for example for human cancer samples.
These assumptions, when valid, permit a much quicker and more powerful significance test to be conducted.

The \code{statistics} are usually a set of probe-wise statistics arising for some comparison from a microarray experiment.
They may be t-statistics, meaning that the genewise null hypotheses would be rejected for large positive or negative values, or they may be F-statistics, meaning that only large values are significant.
Any set of signed statistics, such as log-ratios, M-values or moderated t-statistics, are treated as t-like.
Any set of unsigned statistics, such as F-statistics, posterior probabilities or chi-square tests are treated as F-like.
If \code{type="auto"} then the statistics will be taken to be t-like if they take both positive and negative values and will be taken to be F-like if they are all of the same sign.

There are four possible alternatives to test for.
\code{alternative=="up"} means the genes in the set tend to be up-regulated, with positive t-statistics.
\code{alternative=="down"} means the genes in the set tend to be down-regulated, with negative t-statistics.
\code{alternative=="either"} means the set is either up or down-regulated as a whole.
\code{alternative=="mixed"} test whether the genes in the set tend to be differentially expressed, without regard for direction.
In this case, the test will be significant if the set contains mostly large test statistics, even if some are positive and some are negative.

The latter three alternatives are appropriate if you have a prior expection that all the genes in the set will react in the same direction.
The \code{"mixed"} alternative is appropriate if you know only that the genes are involved in the relevant pathways, possibly in different directions.
The \code{"mixed"} is the only meaningful alternative with F-like statistics.

The test statistic used for the gene-set-test is the mean of the statistics in the set.
If \code{ranks.only} is \code{TRUE} the only the ranks of the statistics are used.
In this case the p-value is obtained from a Wilcoxon test.
If \code{ranks.only} is \code{FALSE}, then the p-value is obtained by simulation using \code{nsim} random sets of genes.
}

\note{
This function does not does allow for inter-gene correlation.
For this reason, the alternative \code{camera} is now recommended over \code{geneSetTest} in those contexts for which \code{camera} is applicable.
}

\seealso{
\code{\link{camera}}, \code{\link{roast}}, \code{\link{romer}}, \code{\link{wilcox.test}}, \code{\link{barcodeplot}}

An overview of tests in limma is given in \link{08.Tests}.
}

\author{Gordon Smyth and Di Wu}

\references{
Goeman, JJ, and Buhlmann P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 

Michaud, J, Simpson, KM, Escher, R, Buchet-Poyau, K, Beissbarth, T, Carmichael, C, Ritchie, ME, Schutz, F, Cannon, P, Liu, M, Shen, X, Ito, Y, Raskind, WH, Horwitz, MS, Osato, M, Turner, DR, Speed, TP, Kavallaris, M, Smyth, GK, and Scott, HS (2008).
Integrative analysis of RUNX1 downstream pathways and target genes.
\emph{BMC Genomics} 9, 363.
\url{http://www.biomedcentral.com/1471-2164/9/363}
}

\examples{
stat <- rnorm(100)
sel <- 1:10; stat[sel] <- stat[sel]+1
wilcoxGST(sel,stat)
}

\keyword{htest}
