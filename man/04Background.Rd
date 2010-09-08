\name{04.Background}
\alias{04.Background}
\title{Background Correction}

\description{
This page deals with background correction methods provided by the \code{\link{backgroundCorrect}}, \code{\link{kooperberg}} or \code{\link{neqc}} functions.
Microarray data is typically background corrected by one of these functions before normalization and other downstream analysis.

\code{backgroundCorrect} works on matrices, \code{EListRaw} or \code{RGList} objects, and calls \code{\link{backgroundCorrect.matrix}}.

The \code{movingmin} method of \code{backgroundCorrect} uses utility functions \code{\link{ma3x3.matrix}} and \code{\link{ma3x3.spottedarray}}.

The \code{normexp} method of \code{backgroundCorrect} uses utility functions \code{\link{normexp.fit}} and \code{\link{normexp.signal}}.

\code{\link{kooperberg}} is a Bayesian background correction tool designed specifically for two-color GenePix data.
It is computationally intensive and requires several additional columns from the GenePix data files.
These can be read in using \code{read.maimages} and specifying the \code{other.columns} argument.

\code{\link{neqc}} is for single-color data.
It performs normexp background correction and quantile normalization using control probes.
It uses utility functions \code{\link{normexp.fit.control}} and \code{\link{normexp.signal}}.
If \code{robust=TRUE}, then \code{normexp.fit.control} uses the function \code{huber} in the MASS package.
}

\author{Gordon Smyth}
\keyword{documentation}
