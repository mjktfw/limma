\title{Multidimensional scaling plot of microarray data}
\name{plotMDS}
\alias{plotMDS}
\description{
Plot the sample relations based on MDS.
}
\usage{
plotMDS(x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1,2), ndim=max(dim.plot),...)
}
\arguments{
  \item{x}{any data object which can be coerced to a matrix, such as \code{ExpressionSet} or \code{EList}.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then defaults the index of the samples.}
  \item{col}{numeric or character vector of colors for the plotting characters.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{ndim}{number of dimensions in which data is to be represented}
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the microarray context is used.
The distance between each pair of samples (columns) is the root-mean-square deviation for the top \code{top} genes which best distinguish that pair of samples.
That is, Euclidean distance is used, but for a different gene subset for each pair of samples.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}

\author{Di Wu and Gordon Smyth}

\seealso{
\code{\link{cmdscale}}

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\examples{
# Simulate gene expression data for 1000 probes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
x <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(x) <- paste("Gene",1:1000)
x[1:50,4:6] <- x[1:50,4:6] + 2
# without labels, indexes of samples are plotted.
plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# with labels as groups, group indicators are plotted.
plotMDS(x,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
}

\keyword{hplot}
