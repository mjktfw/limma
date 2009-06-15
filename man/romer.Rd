\name{romer}
\alias{romer}
\title{Rotation Gene Set Enrichment Analysis}
\description{
Rotation-mean50-rank version of GSEA (gene set enrichment analysis) for linear models.
}
\usage{
romer(iset=NULL,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,nrot=10000)
}
\arguments{
  \item{iset}{list of indices specifying the rows of \code{y} in the gene sets. The list can be made by \link{symbols2indices} and the gene sets can be retrieved from the molecular signatures database in Broad Institute.}
  \item{y}{numeric matrix giving log-expression values.}
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{array.weights}{optional numeric vector of array weights.} 
  \item{block}{optional vector of blocks.}
  \item{correlation}{correlation between blocks.}
  \item{nrot}{number of rotations used to estimate the p-values.}
}
\value{
Matrix with the rows corresponding to estimated p-values for each gene set and the columns corresponding to the number of genes for each gene set and the alternative hypotheses mixed, up, down or either.
}
\details{
This function implements a GSEA of a battery of gene sets similar in motivation to Subramanian et al (2005) but designed for use with linear models.
In is a competitive test, in that the different gene sets are pitted against one another.
The genes are ranked by moderated t-statistic, then the mean rank of the top 50% of genes in each set forms the score for that set.
Instead of permuation, it uses rotation, which is a smoothed version of permutation suitable for linear models (Langsrud, 2005).
\code{romer} can be used with any linear model with some level of replication.

Number of genes is given for each gene set.

p-values are given for four possible alternative hypotheses.
The alternative "up" means the genes in the set tend to be up-regulated, with positive t-statistics.
The alternative "down" means the genes in the set tend to be down-regulated, with negative t-statistics.
The alternative "either" means the set is either up or down-regulated as a whole.
The alternative "mixed" test whether the genes in the set tend to be differentially expressed, without regard for direction.
In this case, the test will be significant if the set contains mostly large test statistics, even if some are positive and some are negative.

The first three alternatives are appropriate if you have a prior expection that all the genes in the set will react in the same direction.
The "mixed" alternative is appropriate if you know only that the genes are involved in the relevant pathways, without knowing the direction of effect for each gene.
The "mixed" alternative is the only one possible with F-like statistics.

Note that \code{romer} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals.
This means that the p-values will vary slightly from run to run.
To get more precise p-values, increase the number of rotations \code{nrot}.
The strategy of random rotations is due to Langsrud (2005).
}

\seealso{
\code{\link{roast}},
\code{\link{geneSetTest}}
}
\author{Yifang Hu and Gordon Smyth}


\references{
Langsrud, O, 2005.
Rotation tests.
\emph{Statistics and Computing} 15, 53-60

Subramanian, A, Tamayo, P, Mootha, VK, Mukherjee, S, Ebert, BL, Gillette, MA, Paulovich, A, Pomeroy, SL, Golub, TR, Lander, ES and Mesirov JP, 2005.
Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
\emph{Proc Natl Acad Sci U S A} 102, 15545-15550
}

\examples{
y <- matrix(rnorm(100*4),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))
iset <- 1:5
y[iset,3:4] <- y[iset,3:4]+3

iset1 <- 1:5
iset2 <- 6:10
romer(iset=list(iset1=iset1,iset2=iset2),y=y,design=design,contrast=2)
}
\keyword{htest}
