\title{modifyWeights}
\name{modifyWeights}
\alias{modifyWeights}
\description{
Modify weights matrix for given gene status values.
}
\usage{
modifyWeights(weights, status, values, multipliers)
}
\arguments{
  \item{weights}{numeric matrix of relative weights, rows corresponding to genes and columns to arrays}
  \item{status}{character vector giving the control status of each spot on the array, of same length as the number of rows of \code{weights}}
  \item{values}{character vector giving subset of the unique values of \code{status}}
  \item{multipliers}{numeric vector of same length as \code{values} giving factor by which weights will be modified}
}

\details{
The function is usually used to temporarily modify the weights matrix during normalization of data.
The function can be used for example to give zero weight to spike-in ratio control spots during normalization.
}

\value{
Numeric matrix of same dimensions as \code{weights} with rows corresponding to \code{values} in \code{status} modified by the specified multipliers.}
\author{Gordon Smyth}
\examples{
MA <- new("MAList")
MA$A <- runif(300,4,16)
MA$M <- rt(300,df=3)
status <- rep("Gene",300)
status[1:3] <- "M=0"
MA$M[1:3] <- 0
status[4:6] <- "M=3"
MA$M[4:6] <- 3
status[7:9] <- "M=-3"
MA$M[7:9] <- -3
plotMA(MA,main="MA-Plot with Simulated Data",status=status,values=c("M=0","M=3","M=-3"),col=c("blue","red","green"))

#  Same as above
attr(status,"values") <- c("M=0","M=3","M=-3")
attr(status,"col") <- c("blue","red","green")
plotMA(MA,main="MA-Plot with Simulated Data",status=status)

#  Same as above
MA$genes$Status <- status
plotMA(MA,main="MA-Plot with Simulated Data")
}
\seealso{
An overview of normalization functions available in LIMMA is given in \link{4.Normalization}.
}
\keyword{hplot}
