\name{alias2Symbol}
\alias{alias2Symbol}
\alias{alias2SymbolTable}
\title{Convert Gene Aliases to Official Gene Symbols}
\description{
Maps gene alias names to official gene symbols.
}
\usage{
alias2Symbol(alias, species = "Hs", expand.symbols = FALSE)
alias2SymbolTable(alias, species = "Hs")
}
\arguments{
  \item{alias}{character vector of gene aliases}
  \item{species}{character string specifying the species.
  Possible values are \code{"Dm"}, \code{"Hs"}, \code{"Mm"} or \code{"Rn"}.}
  \item{expand.symbols}{logical, should those elements of \code{alias} which are already official symbols be expanded if they are aliases for other symbols.}
}
\details{
Aliases are mapped via NCBI Entrez Gene identity numbers using Bioconductor organism packages.
Species are \code{"Dm"} for fly, \code{"Hs"} for human, \code{"Mm"} for mouse and \code{"Rn"} for rat.
The user needs to have the appropriate Bioconductor organism package installed.

\code{alias2Symbol} maps a set of aliases to a set of symbols, without necessarily preserving order.
The output vector may be longer or shorter than the original vector, because some aliases might not be found and some aliases may map to more than one symbol.
\code{alias2SymbolTable} maps each alias to a gene symbol and returns a table with one row for each alias.
If an alias maps to more than one symbol, then the first one found is returned.
}
\value{
Character vector of gene symbols.

\code{alias2SymbolTable} returns a vector of the same length and order as \code{alias}, including \code{NA} values where no gene symbol was found.
\code{alias2Symbol} returns an unordered vector which may be longer or shorter than \code{alias}.
}
\author{Gordon Smyth and Yifang Hu}
\seealso{
This function is often used to assist gene set testing, see
\link{10.GeneSetTests}.
}
\examples{
if(require("org.Hs.eg.db")) alias2Symbol(c("PUMA","NOXA","BIM"))
}

\keyword{gene annotation}

