\name{is.fullrank}
\alias{is.fullrank}

\title{Check for Full Column Rank}

\description{
Test whether a numeric matrix has full column rank.
}

\usage{
is.fullrank(x)
}

\arguments{
\item{x}{a numeric matrix for vector}
}

\value{\code{TRUE} or \code{FALSE}}

\details{
This function is used to check the integrity of design matrices in limma, for example after \link[limma:subsetting]{subsetting} operations.
}

\author{Gordon Smyth}

\examples{
# TRUE
is.fullrank(1)
is.fullrank(cbind(1,0:1))

# FALSE
is.fullrank(0)
is.fullrank(matrix(1,2,2))
}
\keyword{algebra}
