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
  \item{object}{object containing expression intensities from a spotted microarray}
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
\section{Methods}{
A method is defined in limma for \code{object} an object of class
\code{\link[limma:normalizeprintorder-methods]{RGList}}.
There are also default methods.

Type \code{helpMethods("normalizeForPrinter")} at the R prompt for a complete list of classes for which methods are defined.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for normalization is given in \link{4.Normalization}.
}
\keyword{models}
