\name{symbols2indices}
\alias{symbols2indices}
\title{Make Gene Set Indices from Symbols}
\description{
Make a list of gene sets symbols into a list of gene sets indices used to create input for romer function.
}
\usage{
symbols2indices(gmtl.official, symbol)
}
\arguments{
  \item{gmtl.official}{list of gene sets containing official gene symbols. The gene sets can be the ones from the molecular signatures database in Broad Institute.}
  \item{symbol}{vector of symbols for the probesets.}
}
\value{
list of indices.
}
\details{
This function takes a list of official gene sets symbols and returns a list of indices of them according to the vector symbol.
}

\seealso{\code{\link[limma]{romer}}}
\author{Gordon Smyth and Yifang Hu}