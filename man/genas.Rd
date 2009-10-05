\name{genas}
\alias{genas}
\alias{genas}
\title{Genuine Association of Gene Expression Profiles}
\description{
Calculates biological correlation between two gene expression profiles.
}
\usage{
genas(fit, coef=c(1,2))
}
\arguments{
 \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit} and followed by \code{eBayes}}
 \item{coef}{numeric vector of length 2 to indicate which contrasts/columns in the fit object are to be used}
  }
\details{
The biological correlation between the true log fold changes of pairs of comparisons is computed. This method is to be applied when multiple groups (such as treatment groups, mutants or knock-outs) are compared back to the same control group.

This method is an extension of the empirical Bayes method of \code{limma}. It aims to separate the technical correlation, which comes from comparing multiple treatment/mutant/knock-out groups to the same control group, from biological correlation, which is the true correlation of the gene expression profiles between two treatment/mutant/knock-out groups.

The fit object should include only differentially expressed genes. One approach is to calculate the true proportion of differentially expressed genes using \code{convest} on the F p-value produced by \code{lmFit}. Any reasonable set of genes displaying some degree of differential expression should be adequate.
}
\value{
	\code{genas} produces a list with the following components.
	  \item{technical.correlation}{estimate of the technical correlation}
	  \item{covariance.matrix}{estimate of the covariance matrix from which the biological correlation is obtained}
	  \item{biological.correlation}{estimate of the biological correlation}
	  \item{deviance}{the likelihood ratio test statistic used to test whether the biological correlation is equal to 0}
	  \item{p.value}{the p.value associated with \code{deviance}}
 }
 \seealso{
\code{\link{lmFit}}, \code{\link{eBayes}}, \code{\link{contrasts.fit}}
}


\author{Belinda Phipson and Gordon Smyth}
\keyword{regression}
\keyword{models}

\examples{
library(limma)
#  Simulate gene expression data,
#  6 microarrays with 100 genes on each array 
set.seed(2004)
y<-matrix(rnorm(600),ncol=6)

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