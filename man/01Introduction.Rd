\name{01.Introduction}
\docType{package}
\alias{01.Introduction}
\alias{limma}
\alias{limma-package}
\title{Introduction to the LIMMA Package}
\description{
LIMMA is a library for the analysis of gene expression microarray data, especially the use of linear models for analysing designed experiments and the assessment of differential expression.
LIMMA provides the ability to analyse comparisons between many RNA targets simultaneously in arbitrary complicated designed experiments.
Empirical Bayesian methods are used to provide stable results even when the number of arrays is small.
The normalization and data analysis functions are for two-colour spotted microarrays.
The linear model and differential expression functions apply to all microarrays including Affymetrix and other multi-array oligonucleotide experiments.
}

\details{
There are three types of documentation available:
\tabular{ll}{
(1) \tab
The \emph{LIMMA User's Guide} can be reached through the "User
Guides and Package Vignettes" links at the top of the LIMMA
contents page.  The function \code{\link{limmaUsersGuide}} gives
the file location of the User's Guide.\cr
(2) \tab
An overview of limma functions grouped by purpose is contained
in the numbered chapters at the top of the LIMMA contents page,
of which this page is the first.\cr
(3) \tab
The LIMMA contents page gives an
alphabetical index of detailed help topics.\cr
}

The function \code{\link{changeLog}} displays the record of changes to the package.
}

\author{Gordon Smyth}
\references{
Smyth, G. K., Yang, Y.-H., Speed, T. P. (2003). Statistical issues in microarray data analysis. \emph{Methods in Molecular Biology} 224, 111-136.

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume \bold{3}, Article 3.
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}

Smyth, G. K. (2005). Limma: linear models for microarray data.
In: \emph{Bioinformatics and Computational Biology Solutions using R and Bioconductor}.
R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds), Springer, New York, 2005.
}
\keyword{package}
