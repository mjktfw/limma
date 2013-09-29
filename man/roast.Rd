\name{roast}
\alias{roast}
\alias{roast.EList}
\alias{roast.MAList}
\alias{roast.default}
\alias{mroast}
\alias{mroast.EList}
\alias{mroast.MAList}
\alias{mroast.default}
\alias{Roast-class}
\alias{show,Roast-method}
\title{Rotation Gene Set Tests}
\description{
Rotation gene set testing for linear models.
}

\usage{
\method{roast}{EList}(y, index=NULL, design=NULL, contrast=ncol(design), set.statistic="mean",
     gene.weights=NULL, array.weights=NULL, weights=NULL, block=NULL, correlation,
     var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=999)
\method{mroast}{EList}(y, index=NULL, design=NULL, contrast=ncol(design), set.statistic="mean",
     gene.weights=NULL, array.weights=NULL, weights=NULL, block=NULL, correlation,
     var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=999, adjust.method="BH", midp=TRUE, sort="directional")
}

\arguments{
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
  Rows correspond to probes and columns to samples.
  If either \code{var.prior} or \code{df.prior} are null, then \code{y} should contain values for all genes on the arrays. If both prior parameters are given, then only \code{y} values for the test set are required.}
  \item{index}{index vector specifying which rows (probes) of \code{y} are in the test set.  This can be a vector of indices, or a logical vector of the same length as \code{statistics}, or any vector such as \code{y[index,]} contains the values for the gene set to be tested.
  For \code{mroast}, \code{index} is a list of index vectors.}
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{set.statistic}{summary set statistic. Possibilities are \code{"mean"},\code{"floormean"},\code{"mean50"} or \code{"msq"}.}
  \item{gene.weights}{optional numeric vector of weights for genes in the set. Can be positive or negative.  For \code{mroast} this vector must have length equal to \code{nrow(y)}.  For \code{roast}, can be of length \code{nrow(y)} or of length equal to the number of genes in the test set.} 
  \item{array.weights}{optional numeric vector of array weights.}
  \item{weights}{optional matrix of observation weights. If supplied, should be of same dimensions as \code{y} and all values should be positive. If \code{y} is an \code{EList} or \code{MAList} object containing weights, then those weights will be used.}
  \item{block}{optional vector of blocks.}
  \item{correlation}{correlation between blocks.}
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{trend.var}{logical, should a trend be estimated for \code{var.prior}?  See \code{eBayes} for details.  Only used if \code{var.prior} or \code{df.prior} are \code{NULL}.}
  \item{nrot}{number of rotations used to estimate the p-values.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for possible values.}
  \item{midp}{logical, should mid-p-values be used in instead of ordinary p-values when adjusting for multiple testing?}
  \item{sort}{character, whether to sort output table by directional p-values (\code{"directional"}), non-directional p-value (\code{"mixed"}), or not at all (\code{"none"}).}
}

\value{
\code{roast} produces an object of class \code{"Roast"}.
This consists of a list with the following components:
  \item{p.value}{data.frame with columns \code{Active.Prop} and \code{P.Value}, giving the proportion of genes in the set contributing meaningfully to significance and estimated p-values, respectively.
Rows correspond to the alternative hypotheses mixed, up or down.}
  \item{var.prior}{prior value for residual variances.}
  \item{df.prior}{prior degrees of freedom for residual variances.}

\code{mroast} produces a data.frame with a row for each set and the following columns:
  \item{NGenes}{number of genes in set}
	\item{PropDown}{proportion of genes in set with \code{z < -sqrt(2)}}
	\item{PropUp}{proportion of genes in set with \code{z > sqrt(2)}}
	\item{Direction}{direction of change, \code{"Up"} or \code{"Down"}}
	\item{PValue}{two-sided directional p-value}
	\item{FDR}{two-sided directional false discovery rate}
	\item{PValue.Mixed}{non-directional p-value}
	\item{FDR.Mixed}{non-directional false discovery rate}
}

\details{
These functions implement the ROAST gene set tests proposed by Wu et al (2010).
They perform \emph{self-contained} gene set tests in the sense defined by Goeman and Buhlmann (2007).
For \emph{competitive} gene set tests, see \code{\link{camera}}.
For a gene set enrichment analysis style analysis using a database of gene sets, see \code{\link{romer}}.

\code{roast} and \code{mroast} test whether any of the genes in the set are differentially expressed.
They can be used for any microarray experiment which can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model.
If \code{contrast} is not specified, then the last coefficient in the linear model will be tested.

The argument \code{gene.weights} allows directional weights to be set for individual genes in the set.
This is often useful, because it allows each gene to be flagged as to its direction and magnitude of change based on prior experimentation.
A typical use is to make the \code{gene.weights} \code{1} or \code{-1} depending on whether the gene is up or down-regulated in the pathway under consideration.

The arguments \code{array.weights}, \code{block} and \code{correlation} have the same meaning as for the \code{\link{lmFit}} function.
The arguments \code{df.prior} and \code{var.prior} have the same meaning as in the output of the \code{\link{eBayes}} function.
If these arguments are not supplied, they are estimated exactly as is done by \code{eBayes}.

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
By default, \code{mroast} reports ordinary p-values but uses mid-p-values (Routledge, 1994) at the multiple testing stage.
Mid-p-values are probably a good choice when using false discovery rates (\code{adjust.method="BH"}) but not when controlling the family-wise type I error rate (\code{adjust.method="holm"}).
}

\note{
The default setting for the set statistic was changed in limma 3.5.9 (3 June 2010) from \code{"msq"} to \code{"mean"}.
}

\seealso{
\code{\link{camera}}, \code{\link{romer}}, \code{\link{geneSetTest}}.

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
\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments.
\emph{Bioinformatics} 26, 2176-2182.
\url{http://bioinformatics.oxfordjournals.org/content/26/17/2176}
}

\examples{
y <- matrix(rnorm(100*4),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))

# First set of 5 genes contains 3 that are genuinely differentially expressed
index1 <- 1:5
y[index1,3:4] <- y[index1,3:4]+3

# Second set of 5 genes contains none that are DE
index2 <- 6:10

roast(y,index1,design,contrast=2)
mroast(y,index1,design,contrast=2)
mroast(y,list(set1=index1,set2=index2),design,contrast=2)
}
\keyword{htest}
