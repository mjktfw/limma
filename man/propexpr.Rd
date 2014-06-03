\name{propexpr}
\alias{propexpr}
\title{Estimate Proportion of Expressed Probes}
\description{Estimate the proportion of microarray probes which are expressed in each array.}
\usage{
propexpr(x, neg.x=NULL, status=x$genes$Status, labels=c("negative","regular"))
}
\arguments{
  \item{x}{matrix or similar object containing raw intensities for a set of arrays.}
  \item{neg.x}{matrix or similar object containing raw intensities for negative control probes for the same arrays. If \code{NULL}, then negative controls must be provided in \code{x}.}
  \item{status}{character vector giving probe types.}
  \item{labels}{character vector giving probe type identifiers.}
  }
\details{
This function estimates the proportion of expressed in a microarray by utilizing the negative control probes.
Illumina BeadChip arrays contain 750~1600 negative control probes.
The expression profile of these control probes can be saved to a separate file by the Illumina BeadStudio software when using it to output the expression profile for regular probes.
The control probe profile could be re-generated if it was not generated when the regular probe profile was created by BeadStudio.
Other microarray platforms can also use this function to estimate the proportion of expressed probes in each array, provided that they have a set of negative control probes.

\code{labels} can include one or two probe type identifiers. 
Its first element should be the identifier for negative control probes (\code{negative} by default).
If \code{labels} only contains one identifier, then it will be assumed to contain the identifier for negative control probes.
By default, \code{regular} is the identifier for regular probes.
}
\value{
Numeric vector giving the proportions of expressed probes in each array.
}

\references{
Shi, W, de Graaf, C, Kinkel, S, Achtman, A, Baldwin, T, Schofield, L, Scott, H, Hilton, D, Smyth, GK (2010). Estimating the proportion of microarray probes expressed in an RNA sample. Nucleic Acids Research 38, 2168-2176.
}

\author{Wei Shi and Gordon Smyth}

\seealso{ 
Description to the control probes in Illumina BeadChips can be found in \code{\link{read.ilmn}}.
}

\examples{
\dontrun{
x <- read.ilmn(files="sample probe profile.txt",ctrlfiles="control probe profile.txt")
propexpr(x, )
}
}

\keyword{models}