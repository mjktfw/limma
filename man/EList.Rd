\name{EList-class}
\docType{class}
\alias{EList-class}
\alias{EListRaw-class}
\title{Expression List - class}

\description{
A list-based S4 classes for storing expression values (E-values) for a set of one-channel microarrays.
\code{EListRaw} holds expression values on the raw scale.
\code{EList} holds expression values on the log scale, usually after background correction and normalization.
\code{EListRaw} objects are usually created by \code{\link{read.maimages}}.
\code{EList} objects are usually created by \code{\link{normalizeBetweenArrays}}.
}

\section{Components}{
\code{EList} objects can be created by \code{new("EList",E)} where \code{E} is a list.
These classes contains no slots (other than \code{.Data}), but objects should contain a list component \code{E} as follows:
\tabular{ll}{
  \code{E} \tab numeric matrix containing the E-values (raw or log-2 expression ratios).  Rows correspond to probes and columns to arrays.
}
Optional components include:
\tabular{ll}{
  \code{weights} \tab numeric matrix of same dimensions as \code{E} containing relative spot quality weights.  Elements should be non-negative.\cr
  \code{other} \tab list containing other matrices, all of the same dimensions as \code{E}.\cr
  \code{genes} \tab data.frame containing probe information. Should have one row for each probe. May have any number of columns.\cr
  \code{targets} \tab data.frame containing information on the target RNA samples.  Rows correspond to arrays.  May have any number of columns.
}
Valid \code{EList} or \code{EListRaw} objects may contain other optional components, but all probe or array information should be contained in the above components.
}

\section{Methods}{
These classes inherit directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
In addition, \code{EList} objects can be \link[limma:subsetting]{subsetted} and \link[limma:cbind]{combined}.
\code{EList} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{ELists} also inherit a \code{\link[methods]{show}} method from the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}, which means that \code{ELists} will print in a compact way.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} is a more formal class in the Biobase package.
}

\keyword{classes}
\keyword{data}
