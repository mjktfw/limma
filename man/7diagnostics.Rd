\name{7.Diagnostics}
\alias{7.Diagnostics}
\title{Diagnostics and Quality Assessment}

\description{
This page gives an overview of the LIMMA functions available for microarray quality assessment and diagnostic plots.

This package provides an \code{\link[limma:anova-method]{anova}} method which is designed for assessing the quality of an array series or of a normalization method.
It is not designed to assess differential expression of individual genes.
\code{\link[limma:anova-method]{anova}} uses utility functions \code{\link{bwss}} and \code{\link{bwss.matrix}}.

Diagnostics plots can be produced by \code{\link{imageplot}}, \code{\link{plotMA}}, \code{\link{plotPrintTipLoess}}, \code{\link{plotPrintorder}} and \code{\link{plotDensities}}.

\code{plotPrintTipLoess} uses utility functions \code{\link{gridr}} and \code{\link{gridc}}.
\code{\link{plotDensities}} uses utility function \code{\link{RG.MA}}.
}

\author{Gordon Smyth}
\keyword{documentation}
