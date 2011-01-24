\name{roast}
\alias{roast}
\alias{mroast}
\alias{Roast-class}
\alias{show,Roast-method}
\title{Rotation Gene Set Tests}
\description{
Rotation gene set testing for linear models.
}

\usage{
roast(iset=NULL, y, design, contrast=ncol(design), set.statistic="mean",
     gene.weights=NULL, array.weights=NULL, block=NULL, correlation,
     var.prior=NULL, df.prior=NULL, nrot=999)
mroast(iset=NULL, y, design, contrast=ncol(design), set.statistic="mean",
      gene.weights=NULL, array.weights=NULL, block=NULL, correlation,
      var.prior=NULL, df.prior=NULL, nrot=999, adjust.method="BH")
}

\arguments{
  \item{iset}{index vector specifying which rows (probes) of \code{y} are in the test set.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{y[iset,]} contains the values for the gene set to be tested.
  For \code{mroast}, \code{iset} is a list of vectors of indices.}
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
  Rows correspond to probes and columns to samples.
  If either \code{var.prior} or \code{df.prior} are null, then \code{y} should contain values for all genes on the arrays. If both prior parameters are given, then only \code{y} values for the test set are required.}
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{set.statistic}{summary set statistic. Possibilities are \code{"mean"},\code{"floormean"},\code{"mean50"} or \code{"msq"}.}
  \item{gene.weights}{optional numeric vector of weights for genes in the set. Can be positive or negative.} 
  \item{array.weights}{optional numeric vector of array weights.} 
  \item{block}{optional vector of blocks.}
  \item{correlation}{correlation between blocks.}
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{nrot}{number of rotations used to estimate the p-values.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing.
  See \code{\link{p.adjust}} for possible values.}
}

\value{
An object of class \code{"Roast"}.
This consists of a list with the following components:
\tabular{ll}{
  \code{p.value}:\tab
data.frame with columns \code{Active.Prop} and \code{P.Value}, giving the proportion of genes in the set contributing meaningfully to significance and estimated p-values, respectively.
Rows correspond to the alternative hypotheses mixed, up or down.\cr
  \code{var.prior}:\tab prior value for residual variances.\cr
  \code{df.prior}:\tab prior degrees of freedom for residual variances.
}
There is a \code{show} method for \code{"Roast"} objects.
}

\details{
This function implements the ROAST gene set test from Wu et al (2010).
It tests whether any of the genes in the set are differentially expressed.
The function can be used for any microarray experiment which can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model.
If \code{contrast} is not specified, the last coefficient in the linear model will be tested.
The arguments \code{array.weights}, \code{block} and \code{correlation} have the same meaning as they for for the \code{\link{lmFit}} function.

The arguments \code{df.prior} and \code{var.prior} have the same meaning as in the output of the \code{\link{eBayes}} function.
If these arguments are not supplied, they are estimated exactly as is done by \code{eBayes}.

The argument \code{gene.weights} allows directions or weights to be set for individual genes in the set.

The gene set statistics \code{"mean"}, \code{"floormean"}, \code{"mean50"} and \code{msq} are defined by Wu et al (2010).
The different gene set statistics have different sensitivities to small number of genes.
If \code{set.statistic="mean"} then the set will be statistically significantly only when the majority of the genes are differentially expressed.
\code{"floormean"} and \code{"mean50"} will detect as few as 25\% differentially expressed.
\code{"msq"} is sensitive to even smaller proportions of differentially expressed genes, if the effects are reasonably large.

The output gives p-values three possible alternative hypotheses, 
\code{"Up"} to test whether the genes in the set tend to be up-regulated, with positive t-statistics,
\code{"Down"} to test whether the genes in the set tend to be down-regulated, with negative t-statistics,
and \code{"Mixed"} to test whether the genes in the set tend to be differentially expressed, without regard for direction.

\code{roast} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals (Langsrud, 2005), so p-values will vary slightly from run to run.
To get more precise p-values, increase the number of rotations \code{nrot}.
The p-value is computed as \code{(b+1)/(nrot+1)} where \code{b} is the number of rotations giving a more extreme statistic than that observed (Phipson and Smyth, 2010).
This means that the smallest possible p-value is \code{1/(nrot+1)}.

\code{mroast} does roast tests for multiple sets, including adjustment for multiple testing.
}

\note{
The default setting for the set statistic was changed in limma 3.5.9 (3 June 2010) from \code{"msq"} to \code{"mean"}.
}

\seealso{
\code{roast} performs a \emph{self-contained} test in the sense defined by Goeman and Buhlmann (2007).
For a \emph{competitive} gene set test, see \code{\link{wilcoxGST}}.
For a competitive gene set enrichment analysis using a database of gene sets, see \code{\link{romer}.
}

An overview of tests in limma is given in \link{08.Tests}.
}
\author{Gordon Smyth and Di Wu}

\references{
Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 

Langsrud, O (2005).
Rotation tests.
\emph{Statistics and Computing} 15, 53-60.

Phipson B, and Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Article 39.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments.
\emph{Bioinformatics} 26, 2176-2182.
\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq401?}
}

\examples{
y <- matrix(rnorm(100*4),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))

# First set of 5 genes contains 3 that are genuinely differentially expressed
iset1 <- 1:5
y[iset1,3:4] <- y[iset1,3:4]+3

# Second set of 5 genes contains none that are DE
iset2 <- 6:10

roast(iset1,y,design,contrast=2)
mroast(list(set1=iset1,set2=iset2),y,design,contrast=2)
}
\keyword{htest}
