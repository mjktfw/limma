\name{MAList-class}
\docType{class}
\alias{MAList-class}
\title{M-value, A-value Expression List - class}

\description{
A simple list-based class for storing M-values and A-values for a batch of spotted microarrays.
}

\section{Creating Objects from the Class}{
Objects can be created by calls of the form \code{new("MAList",MA)} where \code{MA} is a list.
Objects are normally created by \code{\link{MA.RG}} or \code{\link{normalizeWithinArrays}}.
}

\section{Slots/List Components}{
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

Methods defined in the limma package include \code{\link[limma:normalizebetweenarrays-methods]{normalizeBetweenArrays}},
\code{\link[limma:normalizeprintorder-methods]{normalizeForPrintorder}}
and \code{\link[limma:normalizeWithinArrays-methods]{normalizeWithinArrays}}.

Type \code{showMethods(classes="MAList")} at the R prompt for a complete list of methods for this class.
}

\author{Gordon Smyth}

\seealso{
  \link{2.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[marrayClasses]{marrayNorm-class}} is the corresponding class in the marrayClasses package.
}

\keyword{classes}
\keyword{data}
