\name{4.Normalization}
\alias{4.Normalization}
\title{Normalization of Microarray Data}

\description{
This page gives an overview of the LIMMA functions available to normalize data from spotted two-colour microarrays.
Smyth and Speed (2003) give an overview of the normalization techniques implemented in the functions.

Virtually all data from spotted microarrays will be normalized using \code{\link{normalizeWithinArrays}} with \code{method="printtiploess"}.
A minority of data will also be normalized using \code{\link{normalizeBetweenArrays}} if diagnostic plots suggest a difference in scale between the arrays.

In rare circumstances, data might be normalized using \code{\link{normalizeForPrintorder}} before using \code{\link{normalizeWithinArrays}}.

\code{\link{normalizeWithinArrays}} also does background correction.
At this time, the explicit background correction function \code{\link{backgroundCorrect}} is not used by LIMMA, but the function will become relevant if more sophisticated background correct methods are implemented.
\code{normalizeWithinArrays} uses function \code{\link{normalizeRobustSpline}}.

\code{normalizeBetweenArrays} uses utility functions \code{\link{normalizeMedians}}, \code{\link{normalizeMedianDeviations}} and \code{\link{normalizeQuantiles}}, none of which need to be called directly by users.
}

\author{Gordon Smyth}
\references{
Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. In: \emph{METHODS: Selecting Candidate Genes from DNA Array Screens: Application to Neuroscience}, D. Carter (ed.). To appear. \url{http://www.statsci.org/smyth/pubs/normalize.pdf}
}
\keyword{documentation}
