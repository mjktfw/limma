\name{2.Classes}
\alias{2.Classes}
\title{Classes Defined by this Package}

\description{

This package defines the following data classes.
\describe{
\item{ \code{\link[limma:rglist]{RGList}} }{
	A class used to store raw intensities as they are read in from an image analysis output file,
	usually by \code{read.maimages}.}

\item{ \code{\link[limma:malist]{MAList}} }{
	Intensities converted to M-values and A-values, i.e., to contrasts on the log-scale.
	Usually created from an \code{RGList} using \code{MA.RG} or \code{normalizeBetweenArrays}.}

\item{ \code{\link[limma:marraylm]{MArrayLM}} }{
	Store the result of fitting gene-wise linear models to the normalized intensities or log-ratios.
	Usually created by \code{lmFit}.}
}

Objects of these classes may be \link[limma:subsetting]{subsetted} and multiple data objects may be \link[limma:cbind]{combined}.

All of these classes belong to the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}.
A \code{show} method is defined for \code{LargeDataOject}s which uses the utility function \code{\link{printHead}}.
}

\author{Gordon Smyth}
\keyword{documentation}
