\name{propTrueNull}
\alias{convest}
\alias{propTrueNull}

\title{Estimate Proportion of True Null Hypotheses}

\description{
Estimate the proportion of true null hypotheses from a vector of p-values.
}

\usage{
propTrueNull(p)
convest(p, niter = 100, plot = FALSE, report = FALSE, file="", tol=1e-6)
}

\arguments{
  \item{p}{numeric vector of p-values.}
  \item{niter}{number of iterations to be used in fitting the convex, decreasing density for the p-values.}
  \item{plot}{logical, should updated plots of fitted convex decreasing p-value density be produced at each iteration?}
  \item{report}{logical, should the estimated proportion be printed at each iteration?}
  \item{file}{name of file to which to write the report. Defaults to standard output.}
  \item{tol}{accuracy of the bisectional search for finding a new convex combination of the current iterate and the mixing density}
}

\details{
The proportion of true null hypotheses in a collection of hypothesis tests is often denoted pi0.
\code{propTrueNULL} estimates pi0 from summing local false discovery rates across the p-values.
\code{convest} implements the method of Langaas et al (2005) to estimate pi0 using a convex decreasing density estimate.
}

\value{
Numeric value in the interval [0,1] representing the estimated proportion of true null hypotheses.
}

\references{
Langaas, M., Ferkingstad, E., and Lindqvist, B. (2005).
Estimating the proportion of true null hypotheses, with application
to DNA microarray data.
\emph{Journal of the Royal Statistical Society Series} B, 67, 555-572.
Preprint at \url{http://www.math.ntnu.no/~mettela/pi0.imf}
}

\author{Belinda Phipson and Gordon Smyth for \code{propTrueNull}; Egil Ferkingstad, Mette Langaas and Marcus Davy for \code{convest}}

\seealso{
See \link{08.Tests} for other functions for producing or interpreting p-values.
}

\examples{
# Test statistics
z <- rnorm(200)

# First 40 are have non-zero means
z[1:40] <- z[1:40]+2

# True pi0
160/200

# Two-sided p-values
p <- 2*pnorm(-abs(z))

# Estimate pi0
propTrueNull(p)
convest(p)
}

\keyword{htest}
