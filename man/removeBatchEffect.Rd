\name{removeBatchEffect}
\alias{removeBatchEffect}
\title{Remove Batch Effect}
\description{
Remove batch effects from expression data.
}
\usage{
removeBatchEffect(x, batch=NULL, batch2=NULL, covariates=NULL,
                  design=matrix(1,ncol(x),1), \dots)
}
\arguments{
  \item{x}{numeric matrix, or any data object that can be processed by \code{\link{getEAWP}} containing log-expression values for a series of samples.
  Rows correspond to probes and columns to samples.}
  \item{batch}{factor or vector indicating batches.}
  \item{batch2}{factor or vector indicating batches.}
  \item{covariates}{matrix or vector of covariates to be adjusted for.}
  \item{design}{optional design matrix relating to treatment conditions to be preserved}
  \item{\dots}{other arguments are passed to \code{\link{lmFit}}.}
}
\value{
A numeric matrix of log-expression values with batch and covariate effects removed.
}
\details{
This function is useful for removing batch effects, associated with hybridization time or other technical variables, prior to clustering or unsupervised analysis such as PCA, MDS or heatmaps.
It is not intended to use with linear modelling.
For linear modelling, it is better to include the batch factors in the linear model.

The design matrix is used to describe comparisons between the samples, for example treatment effects, which should not be removed.

The function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects.

The data object \code{x} can be of any class for which \code{lmFit} works.
If \code{x} contains weights, then these will be used in estimating the batch effects.
}

\seealso{
\link{05.Normalization}
}
\author{Gordon Smyth and Carolyn de Graaf}

\examples{
y <- matrix(rnorm(10*6),10,6)
colnames(y) <- c("A1","A2","A3","B1","B2","B3")
y[,1:3] <- y[,1:3] + 10
y
removeBatchEffect(y,batch=c("A","A","A","B","B","B"))
}
