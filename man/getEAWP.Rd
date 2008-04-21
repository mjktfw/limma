\name{getEAWP}
\alias{getEAWP}
\title{Extract Basic Data from Microarray Data Objects}
\description{
Given a microarray data object of any known class, get the expression values, weights, probe annotation and A-values, which are needed for linear modelling.
This function is called by the linear modelling functions in LIMMA.
}
\usage{
getEAWP(object)
}
\arguments{
  \item{object}{a microarray data object.
  An object of class \code{list}, \code{MAList}, \code{marrayNorm}, \code{PLMset}, \code{vsn}, or any class inheriting from \code{ExpressionSet}, or any object that can be coerced to a numeric matrix.}
}
\details{
In the case of two-color objects, the \code{Amean} is computed from the matrix of A-values.
For single-channel objects, \code{Amean} is computed from the matrix of expression vales.
\code{PLMset}, \code{vsn} and \code{ExpressionSet} are assumed to be single-channel for this purpose.

If \code{object} is a matrix, it is assumed to contain log-intensities if the values are all positive and log-ratios otherwise.
\code{Amean} is computed in the former case but not the latter.
}
\value{
A list with components
\item{exprs}{numeric matrix of log-ratios or log-intensities}
\item{weights}{numeric matrix of weights}
\item{probes}{data.frame of probe-annotation}
\item{Amean}{numeric vector of average log-expression for each probe}
}
\author{Gordon Smyth}
\seealso{
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}
