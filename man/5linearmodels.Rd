\name{5.LinearModels}
\alias{5.LinearModels}
\title{Linear Models for Microarrays}

\description{
This page gives an overview of the LIMMA functions available to fit linear models and to interpret the results.

The core of this package is the fitting of gene-wise linear models to microarray data.
The basic idea is to estimate log-ratios between two or more target RNA samples simultaneously.
See the \emph{\link[limma:../doc/usersguide]{LIMMA User's Guide}} for several case studies.
}

\section{Forming the Design Matrix}{
The function \code{\link{designMatrix}} is provided to assist with creation of an appropriate design matrix for two-color microarray experiments using a common reference.
Design matrices for Affymetrix or single-color arrays can be easily created using the ordinary R command \code{\link[base]{model.matrix}}.
For the direct two-color designs the design matrix needs to be created by hand.
}

\section{Fitting Models}{

There are four functions in the package which fit linear models:

\describe{
\item{ \code{\link{lmFit}} }{
	This is a high level function which accepts objects and provides an entry point to the following three functions.}

\item{ \code{\link{lm.series}} }{
	Straightforward least squares fitting of a linear model for each gene.}

\item{ \code{\link{rlm.series}} }{
	An alternative to \code{lm.series} using robust regression as implemented by the \code{rlm} function in the MASS package.}

\item{ \code{\link{gls.series}} }{
	Generalized least squares taking into account correlations between duplicate spots (i.e., replicate spots on the same array).
	The functions \code{\link{duplicateCorrelation}} or \code{\link{dupcor.series}} are used to estimate the inter-duplicate correlation before using \code{gls.series}.}
}
Each of these functions accepts essentially the same argument list and produces a fitted model object of the same form.
The first function \code{lmFit} formally produces an object of class \code{\link[limma:marraylm]{MArrayLM}}.
The other three functions are lower level functions which produce similar output but in unclassed lists.

The main argument is the \bold{design matrix} which specifies which target RNA samples were applied to each channel on each array.
There is considerable freedom to choose the design matrix - there is always more than one choice which is correct provided it is interpreted correctly.
The fitted model object consists of coefficients, standard errors and residual standard errors for each gene.

All the functions which fit linear models use \code{\link{unwrapdups}} which provides an unified method for handling duplicate spots.
}

\section{Making Comparisons of Interest}{

Once a linear model has been fit using an appropriate design matrix, the command \code{\link{makeContrasts}} may be used to form a contrast matrix to make comparisons of interest.
The fit and the contrast matrix are used by \code{\link{contrasts.fit}} to compute fold changes and t-statistics for the contrasts of interest.
This is a way to compute all possible pairwise comparisons between treatments for example in an experiment which compares many treatments to a common reference.
}

\section{Assessing Differential Expression}{

After fitting a linear model, the standard errors are moderated using a simple empirical Bayes model using \code{\link{ebayes}} or \code{\link{eBayes}}.
A moderated t-statistic and a log-odds of differential expression is computed for each contrast for each gene.

\code{\link{ebayes}} and \code{\link{eBayes}} use internal functions \code{\link{fitFDist}}, \code{\link{tmixture.matrix}} and \code{\link{tmixture.vector}}.

The function \code{\link{zscoreT}} is sometimes used for computing z-score equivalents for t-statistics so as to place t-statistics with different degrees of freedom on the same scale.
\code{\link{zcoreGamma}} is used the same way with standard deviations instead of t-statistics.
These functions are for research purposes rather than for routine use.
}

\section{Summarizing Model Fits}{

After the above steps the results may be displayed or further processed using:
\describe{
\item{ \code{\link{toptable}} }{
	Presents a list of the genes most likely to be differentially expressed for a given contrast.}

\item{ \code{\link{classifyTests}} }{
	Uses nested F-tests to classify the genes as up, down or even over the contrasts in the linear model with special attention to genes which are significant in more than one contrast.
	\code{\link{classifyTestsT}} and \code{\link{classifyTestsP}} are simpler methods using cutoffs for the t-statistics or p-values individually.}

\item{ \code{\link{heatdiagram}} }{
	Allows visual comparison of the results across many different conditions in the linear model.
	Not the same as heatdiagrams produced by other packages!
	This function accepts a classification matrix produced by \code{classifyTests}.}

\item{ \code{\link{vennCounts}} }{
	Accepts output from \code{classifyTests} and counts the number of genes in each classification.}

\item{ \code{\link{vennDiagram}} }{
	Accepts output from \code{classifyTests} or \code{vennCounts} and produces a Venn diagram plot.}
}
}

\author{Gordon Smyth}
\references{
Smyth, G. K. (2003). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. \url{http://www.statsci.org/smyth/pubs/ebayes.pdf}

Smyth, G. K., Michaud, J., and Scott, H. (2003). The use of within-array duplicate spots for assessing differential expression in microarray experiments.
\url{http://www.statsci.org/smyth/pubs/dupcor.pdf}
}
\keyword{documentation}
