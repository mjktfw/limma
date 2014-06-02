\name{11.RNAseq}
\alias{11.RNAseq}
\title{Topic: Analysis of RNA-seq Data}

\description{
This page gives an overview of the LIMMA functions available to analyze RNA-seq data.

\describe{
\item{ \code{\link{voom}} }{
	Transform RNA-seq or ChIP-seq counts to log counts per million (log-cpm) with 
associated precision weights.
	After this tranformation, RNA-seq or ChIP-seq data can be analyzed using the same 
functions as would be used for microarray data.}

\item{ \code{\link{diffSplice}} }{
	Test for differential splicing of exons between experimental conditions.}

\item{ \code{\link{topSplice}} }{
	Show a data.frame of top results from \code{diffSplice}.}

\item{ \code{\link{plotSplice}} }{
	Plot results from \code{diffSplice}.}
}
}

\references{
Law, CW, Chen, Y, Shi, W, Smyth, GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\url{http://genomebiology.com/2014/15/2/R29}
}

\keyword{documentation}
