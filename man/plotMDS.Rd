\title{Multidimensional scaling plot of microarray data}
\name{plotMDS}
\alias{plotMDS}
\alias{plotMDS.MDS}
\alias{plotMDS.default}
\alias{MDS-class}
\alias{show,MDS-method}
\description{
Plot the sample relations based on MDS.
Distances on the plot can be interpreted in terms of \emph{leading log2-fold-change}.
}
\usage{
\method{plotMDS}{default}(x, top=500, labels=colnames(x), col=NULL, cex=1, dim.plot=c(1,2), ndim=max(dim.plot), gene.selection="pairwise",
        xlab=paste("Dimension",dim.plot[1]), ylab=paste("Dimension",dim.plot[2]), ...)
\method{plotMDS}{MDS}(x, labels=colnames(x$distance.matrix), col=NULL, cex=1, dim.plot=x$dim.plot, xlab=paste("Dimension",dim.plot[1]), ylab=paste("Dimension",dim.plot[2]),...)

}
\arguments{
  \item{x}{any data object which can be coerced to a matrix, such as \code{ExpressionSet} or \code{EList}.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. If \code{x} has no column names, then defaults the index of the samples.}
  \item{col}{numeric or character vector of colors for the plotting characters.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{which two dimensions should be plotted, numeric vector of length two.}
  \item{ndim}{number of dimensions in which data is to be represented}
  \item{gene.selection}{character, \code{"pairwise"} to choose the top genes separately for each pairwise comparison between the samples or \code{"common"} to select the same genes for all comparisons}
  \item{xlab}{title for the x-axis}
  \item{ylab}{title for the y-axis}
  \item{...}{any other arguments are passed to \code{plot}.}
}

\details{
This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the microarray context is used.
The distance between each pair of samples (columns) is the root-mean-square deviation (Euclidean distance) for the top \code{top} genes.
Distances on the plot can be interpreted as \emph{leading log2-fold-change}, meaning
the typical (root-mean-square) log2-fold-change between the samples for the genes that distinguish those samples.

If \code{gene.selection} is \code{"common"}, then the top genes are those with the largest standard deviations between samples.
If \code{gene.selection} is \code{"pairwise"}, then a different set of top genes is selected for each pair of samples.
The pairwise feature selection may be appropriate for microarray data when different molecular pathways are relevant for distinguishing different pairs of samples.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{
A plot is created on the current graphics device.

An object of class \code{"MDS"} is invisibly returned.
This is a list containing the following components:
\item{distance.matrix}{numeric matrix of pairwise distances between columns of \code{x}}
\item{cmdscale.out}{output from the function \code{cmdscale} given the distance matrix}
\item{dim.plot}{dimensions plotted}
\item{x}{x-xordinates of plotted points}
\item{y}{y-cordinates of plotted points}
\item{gene.selection}{gene selection method}
}

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
mds <- plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# or labels can be provided, here group indicators:
plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
}

\keyword{hplot}
