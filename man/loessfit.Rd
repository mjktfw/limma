\name{loessFit}
\alias{loessFit}
\title{Fast Simple Loess}
\description{
A fast version of \code{loess} when there is only one x-variable and only the fitted values and residuals are required.
}
\usage{
loessFit(y, x, weights=NULL, span=0.3, cell=0.2, iterations=4)
}
\arguments{
  \item{y}{numeric vector of response values.  Missing values are allowed.}
  \item{x}{numeric vector of predictor values  Missing values are allowed.}
  \item{weights}{numeric vector of non-negative weights.  Missing values are allowed.}
  \item{span}{numeric parameter between 0 and 1 specifying proportion of data to be used in loess moving window.
  Larger numbers give smoother fits.}
  \item{cell}{numeric parameter between 0 and 1 used to specify granularity of loess fit}
  \item{iterations}{number of iterations of loess fit}
}

\details{
This function is a low-level equivalent to \code{lowess} in the base library if \code{weights} is null and to \code{loess} in the modreg package otherwise.
It is used by \code{\link{normalizeWithinArrays}}.
The parameters \code{span}, \code{cell} and \code{iterations} have the same meaning as in \code{loess}.
\code{span} is equivalent to the argument \code{f} to \code{lowess} and \code{iterations} is equivalent to \code{iter+1}.
Unlike \code{lowess} this function returns values in original rather than sorted order.

The treatment of missing values is analogous to \code{na.exclude}.
}
\value{
A list with components
\item{fitted}{numeric vector of same length as \code{y} giving the loess fit}
\item{residuals}{numeric vector of same length as \code{x} giving residuals from the fit}
}

\author{Gordon Smyth}

\seealso{
An overview of LIMMA functions for normalization is given in \link{4.Normalization}.

See also \code{\link[base]{lowess}} in the base library and \code{\link[modreg]{loess}} in the modreg package.
}

\examples{
y <- rnorm(1000)
x <- rnorm(1000)
w <- rep(1,1000)
# The following are equivalent apart from execution time
system.time(fit <- loessFit(y,x)$fitted)
system.time(fit <- loessFit(y,x,w)$fitted)
system.time(fit <- fitted(loess(y~x,weights=w,span=0.3,family="symmetric",iterations=4)))
# Similar but with sorted x-values
system.time(fit <- lowess(x,y,f=0.3)$y)
}

\keyword{models}
