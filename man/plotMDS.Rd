\title{Multidimensional scaling plot of microarray data}
\name{plotMDS}
\alias{plotMDS}
\description{
Plot the sample relations based on MDS.
}
\usage{
plotMDS(x,top=500,labels=colnames(x),col=NULL,cex=1,pcs=c(1,2),k=max(pcs),...)
}
\arguments{
  \item{x}{an \code{ExpressionSet}, \code{LumiBatch} object or matrix .}
  \item{top}{number of top genes in used to calculate the Euclidean distance between samples, default to 500.}
  \item{labels}{a character vector describing the different groups or the labels of samples, which will be in the plot. Default is column names of \code{x}. When column names is NA, indexs of samples will be used.}
  \item{col}{a numeric or character vector giving the colors for the plotting characters. Default is NULL in which case default colors in \code{text} will be used. }
  \item{cex}{numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.}
  \item{pcs}{numeric vector of length 2 giving the dimensions of space to visualize the MDS plot, defaults to the first two dimensions.}
  \item{k}{number of dimensions in which data is to be represented}
  \item{...}{any other arguments are passed to \code{plot}}
}

\details{
Estimate the sample relations based on selected genes. Top genes for each distance between two samples are selected by top-Euclidean distance (\code{tEuclid}) measurements. So each distance between two samples has its specific set of top genes in caculation of \code{dist} object.

For any two microarrays, Euclidean distance between the two microarrays for each gene are computed. The top 500 genes with the largest distances are used to caculate the distance between this pair of microarrays. The same procedure is done for each pair of arrays to get the \code{dist} object or a distance matrix. 

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}

\author{Di Wu and Gordon Smyth}

\seealso{
\code{\link{cmdscale}}

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\examples{
# Simulate gene expression data for 1000 probes and 6 microarrays. x is a matrix.
# Microarray are in two groups
# First 50 probes are differentially expressed in second group
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
x <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(x) <- paste("Gene",1:1000)
x[1:50,4:6] <- x[1:50,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# without labels, indexes of samples are plotted.
plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# with labels as groups, group indicators are plotted.
plotMDS(x,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
}

\keyword{hplot}
