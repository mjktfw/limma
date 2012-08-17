\name{loessFit}
\alias{loessFit}
\title{Fast Simple Loess}
\description{
A fast version of locally weighted linear regression when there is only one x-variable and only the fitted values and residuals are required.
}
\usage{
loessFit(y, x, weights=NULL, span=0.3, bin=0.01/(2-is.null(weights)), iterations=4)
}
\arguments{
  \item{y}{numeric vector of response values.  Missing values are allowed.}
  \item{x}{numeric vector of predictor values  Missing values are allowed.}
  \item{weights}{numeric vector of non-negative weights.  Missing values are allowed.}
  \item{span}{numeric parameter between 0 and 1 specifying proportion of data to be used in the local regression moving window.
  Larger numbers give smoother fits.}
  \item{bin}{numeric value between 0 and 1 giving the proportion of the data which can be grouped in a single bin when doing local regression fit.
  \code{bin=0} forces an exact local regression fit with no interpolation.}
  \item{iterations}{number of iterations of loess fit}
}

\details{
This function calls the \code{lowess} function when \code{weights=NULL} and \code{stats:::simpleLoess} otherwise.
Its purpose is to give a unified and streamlined interface to \code{lowess} and \code{loess} for use in \code{\link{normalizeWithinArrays}}.
When \code{weights} is null, this function call \code{lowess} in the stats package with appropropriate choice of tuning parameters.
When \code{weights} is non-null, it is in effect a call to \code{loess} with \code{degree=1}.
See the help pages for those functions for references and credits.

Note that \code{lowess} is faster, needs less memory and is able to use a more accurate interpolation scheme than \code{loess}, so it is desirable to use \code{lowess} whenever \code{loess} is not needed to handle quantitative weights.

The arguments \code{span}, \code{cell} and \code{iterations} here have the same meaning as in \code{loess}.
\code{span} is equivalent to the argument \code{f} of \code{lowess} and \code{iterations} is equivalent to \code{iter+1}.

The parameter \code{bin} is intended to give a uniform interface to the \code{delta} argument of \code{lowess} and the \code{cell} argument of \code{loess}.
\code{bin} translates to \code{delta=bin*diff(range(x))} in a call to \code{lowess} or to \code{cell=bin/span} in a call to \code{loess}.
This is an attempt to put the \code{delta} and \code{cell} arguments on comparable scales.

Unlike \code{lowess}, \code{loessFit} returns values in original rather than sorted order.
Also unlike \code{lowess}, \code{loessFit} allows missing values, the treatment being analogous to \code{na.exclude}.
Unlike \code{loess}, \code{loessFit} returns a linear regression fit if there are insufficient observations to estimate the loess curve.
}
\value{
A list with components
\item{fitted}{numeric vector of same length as \code{y} giving the loess fit}
\item{residuals}{numeric vector of same length as \code{x} giving residuals from the fit}
}

\author{Gordon Smyth}

\seealso{
See \code{\link[stats]{lowess}} and \code{\link[stats]{loess}} in the stats package.

See \link{05.Normalization} for an outline of the limma package normalization functions.
}

\examples{
y <- rnorm(1000)
x <- rnorm(1000)
w <- rep(1,1000)
# The following are equivalent apart from execution time
# and interpolation inaccuracies
fit1 <- loessFit(y,x)$fitted
fit2 <- loessFit(y,x,w)$fitted
fit3 <- fitted(loess(y~x,weights=w,span=0.3,family="symmetric",iterations=4))
# The same but with sorted x-values
fit4 <- lowess(x,y,f=0.3)$y
}

\keyword{models}
