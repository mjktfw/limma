\name{RGList-class}
\docType{class}
\alias{RGList-class}
\title{Red, Green Intensity List - class}

\description{
A simple list-based class for storing red and green channel foreground and background intensities for a batch of spotted microarrays.
\code{RGList} objects are normally created by \code{\link{read.maimages}}.
}

\section{Slots/List Components}{
\code{RGList} objects can be created by \code{new("RGList",RG)} where \code{RG} is a list.
Objects of this class contains no slots (other than \code{.Data}), but objects should contain the following list components:
\tabular{ll}{
  \code{R}:\tab numeric matrix containing the red (cy5) foreground intensities.  Rows correspond to spots and columns to arrays.\cr
  \code{G}:\tab numeric matrix containing the green (cy3) foreground intensities.\cr
  \code{Rb}:\tab numeric matrix containing the red (cy5) background intensities.\cr
  \code{Gb}:\tab numeric matrix containing the green background intensities.\cr
  \code{weights}:\tab optional numeric matric containing relative spot quality weights. Should be non-negative.
}
All of these matrices should have the same dimensions.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.

Methods defined in limma for this class include \code{\link[limma:merge.rglist]{merge}}, \code{\link[limma:normalizebetweenarrays-methods]{normalizeBetweenArrays}},
\code{\link[limma:normalizeprintorder-methods]{normalizeForPrintorder}},
\code{\link[limma:normalizeWithinArrays-methods]{normalizeWithinArrays}}.

\code{RGList} objects inherit a \code{show} method from the virtual class \code{LargeDataObject}.
}

\author{Gordon Smyth}

\seealso{
  \link{2.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[marrayClasses]{marrayRaw-class}} is the corresponding class in the marrayClasses package.
}

\keyword{classes}
\keyword{data}
