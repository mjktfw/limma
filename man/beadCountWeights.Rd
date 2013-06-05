\name{beadCountWeights}
\alias{beadCountWeights}

\title{Bead Count Weights for Illumina BeadChips}

\description{
Estimates weights which account for biological variation and technical variation resulting from varying bead numbers.
}

\usage{
beadCountWeights(y, x, design = NULL, bead.stdev = NULL, nbeads = NULL, array.cv = TRUE, scale = FALSE)
}

\arguments{
 \item{y}{any matrix-like object containing log-expression values that can be coerced to a matrix.}
 \item{x}{any matrix-like object containing raw expression values that can be coerced to a matrix.}
 \item{design}{the design matrix of the microarray experiment, with rows
           corresponding to arrays and columns to coefficients to be
           estimated.  Defaults to the unit vector meaning that the
           arrays are treated as replicates.}
 \item{bead.stdev}{numeric matrix containing bead-level standard deviations.}
 \item{nbeads}{numeric matrix containing number of beads.}
 \item{array.cv}{logical, should technical variation for each observation be calculated from a constant or array-specific coefficient of variation?  The default is to use array-specific coefficients of variation.}
 \item{scale}{logical, should weights be scaled so that the average weight size is the mean of the inverse technical variance along a probe? By default, weights are scaled so that the average weight size along a probe is 1.}
}

\details{
The relative reliability of each gene on each array is measured by estimating its technical and biological variability.

The technical variance for each gene on each array is inversely proportional to the number of beads and is estimated using array-specific bead-level coefficients of variation.

Coefficients of variation are calculated using raw expression values. 

The biological variance for each gene across the arrays are estimated using Newton's iterations, with the assumption that the total residual deviance for each gene from \code{lmFit} is inversely proportional to the sum of the technical variance and biological variance. 

The bead number weights are an inverse of the sum of estimates for technical variances and biological variances.

If any of the arguments \code{design}, \code{bead.stdev} or \code{nbeads} are set explicitly in the call they will over-ride the slots or components in the data \code{object}. The argument \code{design} does not normally need to be set in the call but will be extracted from the data \code{object} if available. If arguments \code{bead.stdev} and \code{nbeads} are not set explicitly in the call, it is necessary that they are available for extraction from the data \code{object}.
}

\value{
A list object with the following components:

\item{weights}{numeric matrix of bead number weights}
\item{cv.constant}{numeric value of constant bead-level coefficient of variation}
\item{cv.array}{numeric vector of array-specific bead-level coefficient of variation}
\item{var.technical}{numeric matrix of technical variance}
\item{var.biological}{numeric vector of biological variance}
}

\references{
Law, CW (2013).
\emph{Precision weights for gene expression analysis}.
PhD Thesis. University of Melbourne, Australia.
}

\author{Charity Law and Gordon Smyth}

\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
