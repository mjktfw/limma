\name{selectModel}
\alias{selectModel}
\title{Select Appropriate Linear Model}
\description{Select the best fitting linear model for each gene by minimizing an information criterion.}
\usage{
selectModel(y, designlist, criterion="bic", var.prior=NULL, df.prior=0, ...)
}
\arguments{
  \item{y}{object of class \code{numeric}, \code{matrix}, \code{MAList}, \code{marrayNorm}, \code{ExpressionSet} or \code{PLMset} containing log-ratios or log-values of expression for a series of microarrays}
  \item{designlist}{list of design matrices} 
  \item{criterion}{information criterion to be used, \code{"bic"}, \code{"aic"} or \code{"lik"}.}
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{...}{other optional arguments to be passed to \code{lmFit}}
}

\value{
List with components
\item{IC}{matrix of information criterion scores, rows for probes and columns for models}
\item{pref}{factor indicating the model with best (lowest) information criterion score}
}

\details{
This function chooses, for each probe, the best fitting model out of a set of models represented by a list of design matrices.
}

\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Alicia Oshlack and Gordon Smyth}
\examples{
nprobes <- 100
narrays <- 5
y <- matrix(rnorm(nprobes*narrays),nprobes,narrays)
A <- c(0,0,1,1,1)
B <- c(0,1,0,1,1)
designlist <- list(
  None=cbind(Int=c(1,1,1,1,1)),
  A=cbind(Int=1,A=A),
  B=cbind(Int=1,B=B),
  Both=cbind(Int=1,AB=A*B),
  Add=cbind(Int=1,A=A,B=B),
  Full=cbind(Int=1,A=A,B=B,AB=A*B)
)
out <- selectModel(y,designlist)
table(out$pref)
}
\keyword{models}
\keyword{regression}
