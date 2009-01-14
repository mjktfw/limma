\name{alias2Symbol}
\alias{alias2Symbol}
\title{Convert Gene Alias to Official Gene Symbols}
\description{
Map a set of gene alias names to a set of gene symbols.
}
\usage{
alias2Symbol(alias, species = "Hs", expand.symbols = FALSE)
}
\arguments{
  \item{alias}{character vector of gene aliases}
  \item{species}{character string specifying the species.
  Possible values are \code{"Dm"}, \code{"Hs"}, \code{"Mm"} or \code{"Rn"}.}
  \item{expand.symbols}{logic, should those aliases which are already official symbols be expanded if they are aliases for other symbols.}
}
\details{
Aliases are mapped via NCBI Entrez Gene identity numbers using Bioconductor organism packages.
}
\value{
A character vector of gene symbols.
May be longer or shorter than the original vector.
}
\author{Gordon Smyth}
\seealso{
\code{\link{unwrapdups}}
}
\examples{
if(!require("org.Hs.eg.db")) alias2Symbol(c("PUMA","NOXA"))
}
\keyword{character}