\name{normalizeForPrintorder}
\alias{normalizeForPrintorder}
\title{Print-Order Normalization - generic}

\description{
Normalize intensity values on one or more spotted microarrays to adjust for print-order effects.
}
\usage{
normalizeForPrintorder(object, ...)
}
\arguments{
  \item{object}{object contain expression intensities from a spotted microarray}
  \item{...}{any other arguments}
}
\details{
Print-order is associated with the 384-well plates used in the printing of spotted microarrays.
There may be variations in DNA concentration or quality between the different plates.
The may be variations in ambient conditions during the time the array is printed.
}
\value{
An object with normalized intensities.
}
\author{Gordon Smyth}
\seealso{
See \code{\showMethods("normalizeForPrinter")} at the R prompt to see what classes of objects have associated methods for this function.
See also \code{\link{normalizeWithinArrays}} and \code{\link{normalizeBetweenArrays}}
}
\keyword{models}
