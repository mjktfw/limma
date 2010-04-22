\name{normalizeCyclicLoess}
\alias{normalizeCyclicLoess}
\title{Normalize Columns of a Matrix by Cyclic Loess}
\description{
Normalize the columns of a matrix, cyclicly applying loess normalization to normalize each pair of columns to each other.
}
\usage{
normalizeCyclicLoess(x, weights = NULL, span=0.4, iterations = 3)
}
\arguments{
  \item{x}{numeric matrix, or object which can be coerced to a numeric matrix, containing log-expression values.}
  \item{weights}{numeric vector of probe weights. Must be non-negative.}
  \item{span}{span of loess smoothing window, between 0 and 1.}
  \item{iterations}{number of times to cycle through all pairs of columns.}
}
\details{
This function is intended to normalize single channel or A-value microarray intensities between arrays.
Cyclic loess normalization is similar effect and intention to quantile normalization, but with some advantages, in particular the ability to incorporate probe weights.

\link[affy]{normalize.loess} in the affy package also does cyclic loess normalization, using probe subsets instead of probe weights.
}
\value{
A matrix of the same dimensions as \code{x} containing the normalized values.
}
\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003). A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.
}
\author{Yunshun (Andy) Chen and Gordon Smyth}
\seealso{
\code{\link{normalizeQuantiles}}.

\link[affy]{normalize.loess} in the affy package.

An overview of LIMMA functions for normalization is given in \link{05.Normalization}.
}  
\keyword{models}
