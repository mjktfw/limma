\name{6.SingleChannel}
\alias{6.SingleChannel}
\title{Single Channel Analysis of Two-Color Microarrays}

\description{
This page gives an overview of the LIMMA functions fit linear models to two-color data in terms of the individual channels.

The function \code{\link{intraspotCorrelation}} estimates the correlation between the two channels.
The regression function code{\link{lmscFit}} takes the correlation as an argument and fits linear models to the two-color data as in terms of the individual log-intensities.
The function \code{\link{array2channel}} converts two-color format target data frames to single channel format, to facilitate the formulation of the design matrix.
}

\author{Gordon Smyth}
\keyword{documentation}
