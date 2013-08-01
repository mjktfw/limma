\name{removeBatchEffect}
\alias{removeBatchEffect}
\title{Remove Batch Effect}
\description{
Remove batch effects from expression data.
}
\usage{
removeBatchEffect(x,batch=NULL,batch2=NULL,covariates=NULL,design=matrix(1,ncol(x),1))
}
\arguments{
  \item{x}{numeric matrix, or any object that can be coerced to a matrix by \code{as.matrix(x)}, containing log-expression intensities for a series of microarrays.
  Rows correspond to probes and columns to arrays. All values must be non-missing.}
  \item{batch}{factor or vector indicating batches.}
  \item{batch2}{factor or vector indicating batches.}
  \item{covariates}{matrix or vector of covariates to be adjusted for.}
  \item{design}{optional design matrix relating to treatment conditions to be preserved}
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
}

\seealso{
\link{05.Normalization}
}
\author{Gordon Smyth and Carolyn de Graaf}

\keyword{regression}
