\name{genas}
\alias{genas}
\alias{genas}
\title{Genuine Association of Gene Expression Profiles}
\description{
Calculates biological correlation between two gene expression profiles.
}
\usage{
genas(fit, coef=c(1,2), subset="all", plot=FALSE, alpha=0.4)
}

\arguments{
 \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit} and followed by \code{eBayes}.}
 \item{coef}{numeric vector of length 2 indicating which columns in the fit object are to be correlated.}
 \item{subset}{character string indicating which subset of genes to use for the correlation analysis.
       Choices are \code{"all"}, \code{"Fpval"}, \code{"p.union"}, \code{"p.int"}, \code{"logFC"} or \code{"predFC"}.}
 \item{plot}{logical, should a plot be produced showing the technical and biological correlation as ellipses?}
 \item{alpha}{numeric value between 0 and 1 determining the transparency of the technical and biological ellipses if a plot is produced.
 \code{alpha=0} indicates fully transparent and \code{alpha=1} indicates fully opague.}
}

\details{
The biological correlation between the log fold changes of pairs of comparisons is computed. This method is to be applied when multiple groups (such as treatment groups, mutants or knock-outs) are compared back to the same control group.

This method is an extension of the empirical Bayes method of \code{limma}. It aims to separate the technical correlation, which comes from comparing multiple treatment/mutant/knock-out groups to the same control group, from biological correlation, which is the true correlation of the gene expression profiles between two treatment/mutant/knock-out groups.

The \code{subset} argument specifies whether and how the fit object should be subsetted. The default is "all", which uses all genes in the fit object to estimate the biological correlation. Only genes that display evidence of differential expression can be used to estimate the biological correlation. The option "Fpval" chooses genes based on how many F p-values are estimated to be truly significant using the method \code{propNotDE}. This should capture genes that display any evidence of differential expression in either of the two contrasts. The options "p.union" and "p.int" are based on the moderated t p-values from both contrasts. From the \code{propNotDE} method an estimate of the number of p-values truly significant in either of the two contrasts can be obtained. "p.union" takes the union of these genes and "p.int" takes the intersection of these genes. The other options, "logFC" and "predFC" subsets on genes that attain a logFC or predFC at least as large as the 90th percentile of the log fold changes or predictive log fold changes on the absolute scale. 

The \code{plot} option is a logical argument that specifies whether or not to plot the log fold changes of the two contrasts. The biological and technical correlations are overlaid on the log fold change versus log fold change scatterplot using transparent ellipses. \code{library(ellipse)} is required to enable the plotting of ellipses. The \code{alpha} argument takes values between 0 and 1 and controls how transparent the ellipses are.
}
\value{
	\code{genas} produces a list with the following components.
	  \item{technical.correlation}{estimate of the technical correlation}
	  \item{biological.correlation}{estimate of the biological correlation}
	  \item{covariance.matrix}{estimate of the covariance matrix from which the biological correlation is obtained}
	  \item{deviance}{the likelihood ratio test statistic used to test whether the biological correlation is equal to 0}
	  \item{p.value}{the p.value associated with \code{deviance}}
	  \item{n}{the number of genes used to estimate the biological correlation} 
 }
 \seealso{
\code{\link{lmFit}}, \code{\link{eBayes}}, \code{\link{contrasts.fit}}
}

\author{Belinda Phipson and Gordon Smyth}

\references{
Phipson, B. (2013).
\emph{Empirical Bayes modelling of expression profiles and their associations}.
PhD Thesis. University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17614}
}

\examples{
library(limma)
#  Simulate gene expression data,
#  6 microarrays with 1000 genes on each array 
set.seed(2004)
y<-matrix(rnorm(6000),ncol=6)

# two experimental groups and one control group with two replicates each
group<-factor(c("A","A","B","B","control","control"))
design<-model.matrix(~0+group)
colnames(design)<-c("A","B","control")

# fit a linear model
fit<-lmFit(y,design)
contrasts<-makeContrasts(A-control,B-control,levels=design)
fit2<-contrasts.fit(fit,contrasts)
fit2<-eBayes(fit2)

# calculate biological correlation between the gene expression profiles of (A vs control) and (B vs control)
genas(fit2)
}