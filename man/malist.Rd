\name{MAList-class}
\docType{class}
\alias{MAList-class}
\title{M-value, A-value Expression List - class}

\description{
A simple list-based class for storing M-values and A-values for a batch of spotted microarrays.
\code{MAList} objects are usually created during normalization by the functions \code{\link{normalizeWithinArrays}} or \code{\link{MA.RG}}.
}

\section{Slots/List Components}{
\code{MAList} objects can be created by \code{new("MAList",MA)} where \code{MA} is a list.
This class contains no slots (other than \code{.Data}), but objects should contain the following list components:
\tabular{ll}{
  \code{M}:\tab numeric matrix containing the M-values or log expression ratios (base 2).  Rows correspond to spots and columns to arrays.\cr
  \code{A}:\tab numeric matrix containing the A-values or average log expression values (base 2).\cr
  \code{weights}:\tab optional numeric matric containing relative spot quality weights.  Should be non-negative.
}
All of the matrices should have the same dimensions.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.

Functions in LIMMA which operate on \code{MAList} objects include
\code{\link{normalizeWithinArrays}},
\code{\link{normalizeBetweenArrays}}
and \code{\link{normalizeForPrintorder}}

\code{MAList} objects inherit a \code{show} method from the virtual class \code{LargeDataObject}.
}

\author{Gordon Smyth}

\seealso{
  \link{2.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[marrayClasses]{marrayNorm-class}} is the corresponding class in the marrayClasses package.
}

\keyword{classes}
\keyword{data}
