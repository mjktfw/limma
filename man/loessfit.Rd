\name{loessFit}
\alias{loessFit}
\title{Univariate Lowess With Prior Weights}

\description{
Univariate locally weighted linear regression allowing for prior weights.
Returns fitted values and residuals.
}

\usage{
loessFit(y, x, weights=NULL, span=0.3, iterations=4, min.weight=1e-5, max.weight=1e5, equal.weights.as.null=TRUE)
}

\arguments{
  \item{y}{numeric vector of response values.  Missing values are allowed.}
  \item{x}{numeric vector of predictor values  Missing values are allowed.}
  \item{weights}{numeric vector of non-negative weights.  Missing values are treated as zero.}
  \item{span}{positive numeric value between 0 and 1 specifying proportion of data to be used in the local regression moving window.
  Larger numbers give smoother fits.}
  \item{iterations}{number of local regression fits. Values greater than 1 produce robust fits.}
  \item{min.weight}{minimum weight. Any lower weights will be reset.}
  \item{max.weight}{maximum weight. Any higher weights will be reset.}
  \item{equal.weights.as.null}{should equal weights be treated as if weights were \code{NULL}, so that \code{lowess} is called? Applies even if all weights are all zero.}
}

\details{
This function is essentially a wrapper function for \code{lowess} and \code{locfit.raw} with added error checking.
The idea is to provide the classic univariate lowess algoritm of Cleveland (1979) but allowing for prior weights and missing values.

The venerable \code{lowess} code is fast, uses little memory and has an accurate interpolation scheme, so it is an advantage to use it when weights are not needed.
This functions calls \code{lowess} when \code{weights=NULL}, but returns values in original rather than sorted order and allows missing values.
The treatment of missing values is analogous to \code{na.exclude}.

When \code{weights} are provided, this function makes repeated calls to \code{locfit.raw} with \code{deg=1}.
Each iteration applies robustifying weights computed from the previous fit, as described by Cleveland (1981).

By default, \code{weights} that are all equal (even all zero) are treated as if they were \code{NULL}.

In principle this function gives analogous results to \code{loess(y~x,weights=weights,span=span,degree=1,family="symmetric")}.
However \code{lowess} and \code{locfit.raw} scale up much more efficiently for long data vectors.

The arguments \code{span} and \code{iterations} here have the same meaning as for \code{loess}.
\code{span} is equivalent to the argument \code{f} of \code{lowess} while \code{iterations} is equivalent to \code{iter+1} for \code{lowess}.
It gives the total number of fits rather than the number of robustifying fits.

When there are insufficient observations to estimate the loess curve, \code{loessFit} returns a linear regression fit.
This mimics the behavior of \code{lowess} but not that of \code{loess} or \code{locfit.raw}.
}

\value{
A list with components
\item{fitted}{numeric vector of same length as \code{y} giving the loess fit}
\item{residuals}{numeric vector of same length as \code{x} giving residuals from the fit}
}

\author{Gordon Smyth}

\references{
Cleveland, W. S. (1979).
Robust locally weighted regression and smoothing scatterplots.
\emph{Journal of the American Statistical Association} 74, 829-836.
}

\seealso{
This function uses \code{\link{lowess}} and \code{\link[locfit]{locfit.raw}}.
See the help pages of those functions for references and credits.

Compare with \code{\link{loess}} in the stats package.

See \link{05.Normalization} for an outline of the limma package normalization functions.
}

\examples{
x <- (1:100)/101
y <- sin(2*pi*x)+rnorm(100,sd=0.4)
out <- loessFit(y,x)
plot(x,y)
lines(x,out$fitted,col="red")

# Example using weights

y <- x-0.5
w <- rep(c(0,1),50)
y[w==0] <- rnorm(50,sd=0.1)
pch <- ifelse(w>0,16,1)
plot(x,y,pch=pch)
out <- loessFit(y,x,weights=w)
lines(x,out$fitted,col="red")
}

\keyword{models}
