\name{normalizeQuantiles}
\alias{normalizeQuantiles}
\title{Normalize Columns of a Matrix to have the same Quantiles}
\description{
Normalize the columns of a matrix to have the same quantiles, allowing for missing values.
Users do not normally need to call this function directly - use \code{\link{normalizeBetweenArrays}} instead.
}
\usage{
normalizeQuantiles(A, ties=FALSE)
}
\arguments{
  \item{A}{numeric matrix. Missing values are allowed.}
  \item{ties}{logical. If \code{TRUE}, ties in each column of \code{A} are treated in careful way. tied values will be normalized to the mean of the corresponding pooled quantiles.}
}
\details{
This function is intended to normalize single channel or A-value microarray intensities between arrays.
Each quantile of each column is set to the mean of that quantile across arrays.
If \code{ties=TRUE}, then the quantiles are further averaged across the quantiles corresponding to any set of tied values in any column.
If there are no missing values and no ties then the normalized columns are permutations of one another.
}
\value{
A matrix of the same dimensions as \code{A} containing the normalized values.
}
\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003), A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for normalization is given in \link{4.Normalization}.
}  
\keyword{models}
