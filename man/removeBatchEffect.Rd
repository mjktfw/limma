\name{removeBatchEffect}
\alias{removeBatchEffect}
\title{Remove Batch Effect}
\description{
Remove a batch effect from expression data.
}
\usage{
removeBatchEffect(x,batch,design=NULL)
}
\arguments{
  \item{x}{numeric matrix containing log-expression intensity for a series of microarrays. Each column corresponds to an array.}
  \item{batch}{a factor or vector indicating batches.}
  \item{design}{an optional design matrix}
}
\value{
A numeric matrix of log-expression values with batch effects removed.
}
\details{
This function is useful for removing batch effects, associated with hybridization time or other technical variables, prior to clustering or unsupervised analysis.

This function is not intended to use prior to linear modelling.
If linear modelling is intended, it is better to include the batch effect as part of the linear model.
}

\seealso{
\link{05.Normalization}
}
\author{Gordon Smyth and Carolyn de Graaf}

\keyword{regression}
