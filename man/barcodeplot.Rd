\name{barcodeplot}
\alias{barcodeplot}
\title{Barcode Plot}
\description{
Plot positions of one or two gene sets in a ranked list of statistics.
}
\usage{
barcodeplot(statistics, index, index2=NULL, labels=c("Up","Down"), quantiles=c(-1,1), col.bars=NULL, worm.window.width=0, ...)
}
\arguments{
  \item{index}{index vector for the gene set.  This can be a vector of indices, or a logical vector of the same length as \code{statistics} or, in general, any vector such that \code{statistic[index]} gives the statistic values for the gene set to be tested.}
  \item{index2}{index vector for a second gene set.  Usually used to specify down-regulated genes when \code{index} is used for up-regulated genes.vector specifying the elements of \code{statistic} in the test group.}
  \item{statistics}{numeric vector giving the values of statistics to rank genes by.}
  \item{labels}{character vector of labels for the two groups of RNA samples that are being compared by the statistics.  First label is associated with high statistics and is displayed at the left end of the plot.  Second label is associated with low or negative statistics and is displayed at the right end of the plot.}
  \item{quantiles}{numeric vector of length 2, giving cutoff values for \code{statistics} considered small or large respectively.  Used to color the rectangle of the barcodeplot.}
  \item{col.bars}{character vector giving colors for the bars on the barcodeplot. Defaults to \code{"black"} for one set or \code{c("red","blue")} for two sets.}
  \item{worm.window.width}{integer giving each sliding window width for the enrichment worm. The worm shows the relative enrichment of the gene set in each window. The default of \code{0} corresponds to no worm.}
  \item{\ldots}{other arguments are passed to \code{plot}.}
}

\value{
No value is returned but a plot is produced as a side effect.
}

\details{
Rank statistics left to right from largest to smallest, and show positions of one or two specified subsets.
This plot is typically used in conjunction with gene set tests.
It first appeared in the literature in Lim et al (2009).
It was inspired by the set location plot of Subramanian et al (2005).
The enrichment score is calculated by \code{(sliding.window.sum/worm.window.width)/(gene.set.size/statistics.size)}.
}

\seealso{
\code{\link{camera}}, \code{\link{wilcox.test}}
}

\author{Gordon Smyth, Di Wu and Yifang Hu}

\references{
Lim E, Vaillant F, Wu D, Forrest NC, Pal B, Hart AH, Asselin-Labat ML, Gyorki DE, Ward T, Partanen A, Feleppa F, Huschtscha LI, Thorne HJ; kConFab; Fox SB, Yan M, French JD, Brown MA, Smyth GK, Visvader JE, Lindeman GJ (2009).
Aberrant luminal progenitors as the candidate target population for basal tumor development in BRCA1 mutation carriers.
\emph{Nat Med}, 15, 907-913..

Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP (2005).
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
\emph{Proc Natl Acad Sci USA}, 102, 15545-15550.
}

\examples{
stat <- rnorm(100)
sel <- 1:10
stat[sel] <- stat[sel]+1
barcodeplot(stat,sel)
sel2 <- 11:20
stat[sel2] <- stat[sel2]-1
barcodeplot(stat,sel,sel2)
barcodeplot(stat,sel,sel2,worm.window.width = 20)
}
