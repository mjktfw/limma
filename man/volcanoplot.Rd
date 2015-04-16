\title{Volcano Plot}
\name{volcanoplot}
\alias{volcanoplot}
\description{
Creates a volcano plot of log-fold changes versus log-odds of differential expression.
}
\usage{
volcanoplot(fit, coef=1, highlight=0, names=fit$genes$ID,
            xlab="Log Fold Change", ylab="Log Odds", pch=16, cex=0.35, \dots)
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted linear model object}
  \item{coef}{integer giving the coefficient}
  \item{highlight}{number of top genes to be highlighted}
  \item{names}{character vector giving text labels for the probes to be used in highlighting}
  \item{xlab}{character string giving label for x-axis}
  \item{ylab}{character string giving label for y-axis}
  \item{pch}{vector or list of plotting characters. Default is integer code 16 which gives a solid circle.}
  \item{cex}{numeric vector of plot symbol expansions. Default is 0.35.}
  \item{\dots}{any other arguments are passed to \code{plot}}
}

\details{
A volcano plot is any plot which displays fold changes versus a measure of statistical significance of the change.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of presentation plots following the fitting of a linear model in LIMMA is given in \link{06.LinearModels}.
}
\examples{
#  See lmFit examples
}
