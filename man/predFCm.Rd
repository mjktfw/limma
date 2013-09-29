\name{predFCm}
\alias{predFCm}
\alias{predFCm}
\title{Predictive log fold change for microarrays}
\description{
Calculates the predictive log fold change for a particular coefficient from a fit object.
}
\usage{
predFCm(fit,coef=2,prob=TRUE,VarRel=NULL)
}
\arguments{
 \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit} and followed by \code{eBayes}}
 \item{coef}{integer indicating which contrasts/columns in the fit object are to be used}
 \item{prob}{logical, whether the probability that the gene is differentially expressed should be taken into account}
 \item{VarRel}{string, options are "Independent" or "Increasing"}
  }
\details{
The predictive log fold changes are calculated as the posterior log fold changes in the empirical Bayes hierarchical model. The log fold changes are shrunk towards zero depending on how variable they are. The \code{VarRel} argument specifies whether the prior belief is that the log fold changes are independent of the variability of the genes (option "Independent"), or whether the log fold changes increase with increasing variability of the genes (option "Increasing"). The \code{prob} argument is a logical argument indicating whether the probability that a particular gene is differentially expressed should be taken into account.
}
\value{
	\code{predFCm} produces a numeric vector corresponding to the predictive or posterior log fold changes of the specified contrast
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

# output predictive log fold changes for first contrast for first 5 genes
pfc<-predFCm(fit2,coef=1,prob=FALSE)
cbind(pfc[1:5],fit2$coeff[1:5,1])
}