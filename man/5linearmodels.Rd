\name{5.LinearModels}
\alias{5.LinearModels}
\title{Linear Models for Microarrays}

\description{
This page gives an overview of the LIMMA functions available to fit linear models and to interpret the results.

The core of this package is the fitting of gene-wise linear models to microarray data.
The basic idea is to estimate log-ratios between two or more target RNA samples simultaneously.
See the \emph{\link[limma:../doc/usersguide]{LIMMA User's Guide}} for several case studies.
There are three functions in the package which fit linear models:

\describe{
\item{ \code{\link{lm.series}} }{
	Straightforward least squares fitting of a linear model for each gene.}

\item{ \code{\link{rlm.series}} }{
	An alternative to \code{lm.series} using robust regression as implemented by the \code{rlm} function in the MASS package.}

\item{ \code{\link{gls.series}} }{
	Generalized least squares taking into account correlations between duplicate spots (i.e., replicate spots on the same array).
	The function \code{\link{dupcor.series}} is used to estimate the inter-duplicate correlation before using \code{gls.series}.}
}
Each of these functions accepts essentially the same argument list and produces a fitted model object of the same form.
The main argument is the \bold{design matrix} which specifies the comparisons between the RNA targets which you want to estime.
The fitted model object consists of coefficients, standard errors and residual standard errors for each gene.

If the design matrix has orthogonal columns, then the output object can be explanded using \code{\link{contrasts.fit}} to include contrasts not explicitly in the design matrix.
This is a way to compute all possible pairwise comparisons between treatments for example in an experiment which compares many treatments to a common reference.

After fitting a linear model, the standard errors are moderated using a simple empirical Bayes model using \code{\link{ebayes}}.
A log-odds of differential expression is computed for each contrast for each gene.

After linear modelling and moderation, the results may be displayed or further processed using:
\describe{
\item{ \code{\link{toptable}} }{
	Presents a list of the genes most likely to be differentially expressed for a given contrast.}

\item{ \code{\link{heatdiagram}} }{
	Allows visual comparison of the results across the contrasts in the linear model.}

\item{ \code{\link{classifyTests}} }{
	Classify the genes as up, down or even over the contrasts in the linear model.}
}
}

\section{Utility Functions}{
All the functions which fit linear models use \code{\link{unwrapdups}} which gives an unified method for treating duplicate spots.

\code{\link{ebayes}} uses internal functions \code{\link{fitFDist}}, \code{\link{tmixture.matrix}} and \code{\link{tmixture.vector}}.

\code{\link{classifyTests}} uses internal function \code{\link{classifyTests.vector}}.
}
\author{Gordon Smyth}
\keyword{documentation}
