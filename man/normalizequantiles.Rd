\name{normalizeQuantiles}
\alias{normalizeQuantiles}
\title{Normalize Columns of a Matrix to have the same Quantiles}
\description{
Normalize the columns of a matrix to have the same quantiles, allowing for missing values.
Users do not normally need to call this function directly - use \code{normalizeBetweenArrays} instead.
}
\usage{
normalizeQuantiles(A)
}
\arguments{
  \item{A}{numeric matrix. Missing values are allowed.}
}
\details{
This function is intended to normalize single channel or A-value microarray intensities between arrays.
Each quantile of each column is set to the mean of that quantile across arrays.
If there are no missing values then the normalized columns are permutations of one another.
}
\value{
A matrix of the same dimensions as \code{A} containing the normalized values.
}
\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003), A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.
}
\author{Gordon Smyth}
\seealso{
  \code{\link{normalizeBetweenArrays}}, \code{\link{normalizeMedianDeviations}}
}  
\keyword{models}
