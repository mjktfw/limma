\name{removeBatchEffect}
\alias{removeBatchEffect}
\title{Remove Batch Effect}
\description{
Remove a batch effect from expression data.
}
\usage{
removeBatchEffect(x,batch,design=matrix(1,ncol(x),1))
}
\arguments{
  \item{x}{numeric matrix containing log-expression intensity for a series of microarrays. Each column corresponds to an array.}
  \item{batch}{a factor or vector indicating batches.}
  \item{design}{design matrix}
}
\value{
A numeric matrix of log-expression values with batch effects removed.
}
\details{
This function is useful for removing batch effects, associated with hybridization time or other technical variables, prior to clustering or unsupervised analysis.

This function is intended for use with clustering or PCA, not for use prior to linear modelling.
If linear modelling is intended, it is better to include the batch effect as part of the linear model.

The design matrix is used to describe comparisons between the samples, for example treatment effects, which should not be removed.
The batch effect is estimated after adjusting for the design factors.
By default, this is just the grand mean.
}

\seealso{
\link{05.Normalization}
}
\author{Gordon Smyth and Carolyn de Graaf}

\keyword{regression}
