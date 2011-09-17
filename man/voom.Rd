\name{voom}
\alias{voom}
\title{Variance Observational-Level Modelling}
\description{
Transform count data to log2-counts per million, estimate the mean-variance relationship and use this to compute appropriate observational-level weights.
The data are then ready for linear modelling.
}

\usage{
voom(counts, design = NULL, lib.size = NULL, normalize.method = "none", plot = FALSE, ...)
}
\arguments{
 \item{counts}{object of class \code{matrix} containing raw counts.}
 \item{design}{design matrix with rows corresponding to samples and columns to coefficients to be estimated.  Defaults to the unit vector meaning that samples are treated as replicates.}
 \item{lib.size}{numeric vector containing total librazy sizes for each sample. If \code{NULL}, library sizes are calculated as column sums of \code{counts}.}
 \item{normalize.method}{normalization method to be applied to the log2-counts-per-million. Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.}
 \item{plot}{\code{logical} value indicating whether a plot of mean-variance trend is displayed.}
 \item{...}{other arguments are passed to \code{lmFit}.}
  }

\details{
Typically, mean-variance relationships exist in count data. Weights can be assigned to each observation in a linear modelling process to adjust for heteroscedasticity. 

In an experiment, a count value is observed for each tag in each sample. A tag-wise mean-variance trend is computed using \code{link{lowess}}. The tag-wise mean is the mean log2 count with an offset of 0.5, across samples for a given tag. The tag-wise variance is the quarter-root-variance of normalized log2 counts per million values with an offset of 0.5, across samples for a given tag. Tags with zero counts across all samples are not included in the lowess fit. Normalization is performed using \code{\link{normalizeBetweenArrays}}. 

Using fitted values of log2 counts from a linear model fit by \code{\link{lmFit}}, variances from the mean-variance trend were interpolated for each observation. This was carried out by \code{\link{approxfun}}. Inverse variance weights can be used to correct for mean-variance trend in the count data.
}
\value{
An EList object with the following components:

\item{E}{numeric matrix of normalized expression values on the log2 scale}
\item{weights}{numeric matrix of inverse variance weights}
\item{design}{numeric matrix of experimental design}
\item{lib.size}{numeric vector of total library sizes}
 }
\seealso{
\code{\link{normalizeBetweenArrays}}
}

\author{Charity Law and Gordon Smyth}
\keyword{weights}
\keyword{variance}
