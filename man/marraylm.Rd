\name{MArrayLM-class}
\docType{class}
\alias{MArrayLM-class}
\title{Microarray Linear Model Fit - class}

\description{
A class for storing the results of fitting gene-wise linear models to a batch of microarrays.
Objects are normally created by \code{\link{lmFit}}.
}

\section{Slots}{
  \describe{
    \item{\code{genes}:}{\code{data.frame} containing gene names and annotation}
    \item{\code{design}:}{design \code{matrix}}
    \item{\code{contrasts}:}{\code{matrix} defining contrasts of coefficients for which results are desired}
    \item{\code{coefficients}:}{\code{matrix} containing fitted coefficients}
    \item{\code{stdev.unscaled}:}{\code{matrix} containing unscaled standard deviations of the coefficients}
    \item{\code{s2.residual}:}{\code{numeric} vector containing residual variances for each gene}
    \item{\code{df.residual}:}{\code{numeric} vector containing residual degrees of freedom for each gene}
    \item{\code{s2.prior}:}{\code{numeric} value giving empirical Bayes estimated prior value for residual variances}
    \item{\code{df.prior}:}{\code{numeric} vector giving empirical Bayes estimated degrees of freedom associated with \code{s2.prior} for each gene}
    \item{\code{s2.post}:}{\code{numeric} vector giving posterior residual variances}
    \item{\code{tstat}:}{\code{matrix} containing empirical Bayes t-statistics}
    \item{\code{varcoef.prior}:}{\code{numeric} vector giving empirical Bayes estimated variance for each true coefficient}
  }
}

\section{Methods}{
\code{MAList} objects inherit a \code{show} method from the virtual class \code{LargeDataObject}.
}

\author{Gordon Smyth}

\seealso{
  \link{2.Classes} gives an overview of all the classes defined by this package.
}

\keyword{classes}
\keyword{regression}
