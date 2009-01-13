\name{roast}
\alias{roast}
\title{roast}
\description{
Rotation gene set testing for linear models.
}
\usage{
roast(iset=NULL,y,design,contrast=ncol(design),gene.weights=NULL,array.weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,nrot=1000)
}
\arguments{
  \item{iset}{vector specifying the rows of \code{y} in the test set.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{y[selected,]} contains the values for the gene set to be tested.}
  \item{y}{numeric matrix giving log-expression values. If \code{var.prior} or \code{df.prior} are null, then \code{y} should contain values for all genes on the arrays. If both prior parameters are given, then only \code{y} values for the test set are required.}
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{gene.weights}{optional numeric vector of weights for genes in the set.} 
  \item{array.weights}{optional numeric vector of array weights.} 
  \item{block}{optional vector of blocks.}
  \item{correlation}{correlation between blocks.}
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{nrot}{number of rotations used to estimate the p-values.}
}
\value{
data.frame with columns \code{Z}, \code{Active} and \code{P.Value}.
The \code{Z} column gives average (root mean square) z-statistics for the genes in the set.
The \code{Active} gives the proportion of genes in the set contributing meaningfully to significance, defined as those with squared z-values greater than 2.
The \code{P.Value} gives estimated p-values.
The rows correspond to the alternative hypotheses mixed, up, down or either.
}
\details{
This function tests whether any of the genes in the set are differentially expressed.
It uses rotation, which is a smoothed version of permutation suitable for linear models (Langsrud, 2005).
It can be used for any linear model with replication.
and negative values and otherwise will be taken to be F-like.

This is a self-contained test is the sense that genes outside the test set do not play a role
(Goeman, JJ, and Buhlmann P, 2007).
A competitive gene set test is performed by \code{\link{geneSetTest}}.

p-values are given for four possible alternative hypotheses.
\code{alternative=="up"} means the genes in the set tend to be up-regulated, with positive t-statistics.
\code{alternative=="down"} means the genes in the set tend to be down-regulated, with negative t-statistics.
\code{alternative=="either"} means the set is either up or down-regulated as a whole.
\code{alternative=="mixed"} test whether the genes in the set tend to be differentially expressed, without regard for direction.
In this case, the test will be significant if the set contains mostly large test statistics, even if some are positive and some are negative.

The first three alternatives are appropriate if you have a prior expection that all the genes in the set will react in the same direction.
The \code{"mixed"} alternative is appropriate if you know only that the genes are involved in the relevant pathways, without knowing the direction of effect for each gene.
The \code{"mixed"} alternative is the only one possible with F-like statistics.

Note that \code{roast} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals.
This means that the p-values will vary slightly from run to run.
To get more precise p-values, increase the number of rotations \code{nrot}.
The strategy of random rotations is due to Langsrud (2005).
}

\seealso{\code{\link[limma]{geneSetTest}}}
\author{Gordon Smyth and Di Wu}

\references{
Goeman, JJ, and Buhlmann P, 2007.
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 

Langsrud, O. (2005).
Rotation tests.
\emph{Statistics and Computing} 15, 53-60
}

\examples{
y <- matrix(rnorm(100*4),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))
iset <- 1:5
y[iset,3:4] <- y[iset,3:4]+3
roast(iset,y,design,contrast=2)

# Alternative approach useful if multiple gene sets are tested:
fit <- lmFit(y,design)
sv <- squeezeVar(fit$sigma^2,df=fit$df.residual)
iset1 <- 1:5
iset2 <- 6:10
roast(y=y[iset1,],design=design,contrast=2,var.prior=sv$var.prior,df.prior=sv$var.prior)
roast(y=y[iset2,],design=design,contrast=2,var.prior=sv$var.prior,df.prior=sv$var.prior)
}
\keyword{htest}
