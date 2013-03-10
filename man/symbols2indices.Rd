\name{symbols2indices}
\alias{symbols2indices}
\title{Convert Gene Set Symbols to Indices}
\description{
Make a list of gene set symbols into a list of gene sets indices. 
}
\usage{
symbols2indices(gene.sets, symbols, remove.empty=TRUE)
}
\arguments{
  \item{gene.sets}{list of character vectors, each vector containing the symbols for a set of genes.}
  \item{symbols}{character vector of gene symbols.}
  \item{remove.empty}{logical, should sets of size zero be removed from the output?}
}
\value{
list of integer vectors, each vector containing the indices of a gene set in the vector \code{symbols}.
}
\details{
This function used to create input for \code{romer} function.
Typically, \code{symbols} is the vector of symbols of genes on a microarray, and \code{gene.sets} is obtained constructed from a database of gene sets, for example a representation of the Molecular Signatures Database (MSigDB) downloaded from \url{http://bioinf.wehi.edu.au/software/MSigDB}.
}

\seealso{
\code{\link{romer}}, \code{\link{mroast}}, \code{\link{camera}}
}
\author{Gordon Smyth and Yifang Hu}