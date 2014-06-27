\title{Plot of regularized linear discriminant functions for microarray data}
\name{plotRLDF}
\alias{plotRLDF}
\description{
Plot of regularized linear discriminant functions for microarray data.
}
\usage{
plotRLDF(y,design=NULL,z=NULL,labels.y=NULL,labels.z=NULL,col.y=1,col.z=1,
df.prior=5,show.dimensions=c(1,2),main=NULL,nprobes=500,\dots)}
\arguments{
 \item{y}{any data object which can be coerced to a matrix, such as \code{ExpressionSet} or \code{EList}. The training dataset.}
 \item{z}{any data object which can be coerced to a matrix, such as \code{ExpressionSet} or \code{EList}. The dataset to be classified.}
 \item{design}{the design matrix ofthe microarray experiment for \code{y}, with rows corresponding to arrays and columns to coefficients to be estimated. Defaults to the unit vector meaning that the arrays are treated as replicates.}
 \item{labels.y}{character vector of sample names or labels in \code{y}. Default is integers starting from 1.}
 \item{labels.z}{character vector of sample names or labels in \code{z}. Default is \code{letters}.}
  \item{col.y}{numeric or character vector of colors for the plotting characters of \code{y}. Default is black.}
 \item{col.z}{numeric or character vector of colors for the plotting characters of \code{z}. Default is black.}
 \item{df.prior}{prior degrees of freedom for residual variances. Used in gene selection.}
 \item{show.dimensions}{which two dimensions should be plotted, numeric vector of length two.}
 \item{main}{title of the plot.}
 \item{nprobes}{number of probes to be used for the calculations. Selected by moderated F tests.}
 \item{\dots}{any other arguments are passed to \code{plot}.}
}
\details{
This function is a variation on the plot of usual linear discriminant fuction, in that the within-group covariance matrix is regularized to ensure that it is invertible, with eigenvalues bounded away from zero.
A diagonal regulation using \code{df.prior} and the median within-group variance is used.

The calculations are based on a filtered list of probes.
The \code{nprobes} probes with largest moderated F statistics are used to discriminate.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}
\value{A list containing metagene information is (invisibly) returned.
A plot is created on the current graphics device.}

\author{Di Wu and Gordon Smyth}

\seealso{
\code{lda} in package \code{MASS}
}

\examples{

# Simulate gene expression data for 1000 probes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
y <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(y) <- paste("Gene",1:1000)
y[1:50,4:6] <- y[1:50,4:6] + 2

z <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(z) <- paste("Gene",1:1000)
z[1:50,4:6] <- z[1:50,4:6] + 1.8
z[1:50,1:3] <- z[1:50,1:3] - 0.2

design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digit=3)

plotRLDF(y,z, design=design)
}

\keyword{hplot}

