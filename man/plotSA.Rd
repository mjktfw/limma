\title{Sigma vs A plot for microarray linear model}
\name{plotSA}
\alias{plotSA}
\description{
Plot log residual standard deviation versus average log expression for a fitted microarray linear model.
}
\usage{
plotSA(fit, xlab="Average log-expression", ylab="log2(sigma)",
       zero.weights=FALSE, pch=16, cex=0.2, \dots)
}
\arguments{
  \item{fit}{an \code{MArrayLM} object.}
  \item{xlab}{character string giving label for x-axis}
  \item{ylab}{character string giving label for y-axis}
  \item{pch}{vector or list of plotting characters. Default is integer code 16 which gives a solid circle.}
  \item{cex}{numeric expansion factor for plotting character. 
  Defaults to 0.2.}
  \item{zero.weights}{logical, should spots with zero or negative weights be plotted?}
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
This plot is used to check the mean-variance relationship of the expression data, after fitting a linear model.

See \code{\link[graphics]{points}} for possible values for \code{pch} and \code{cex}.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}
