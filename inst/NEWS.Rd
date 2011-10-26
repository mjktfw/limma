\name{limmanews}
\title{CQN News}
\encoding{UTF-8}

\section{Version 3.10.0}{
  \itemize{
    \item New function voom() provides functionality to analyse RNA-Seq data. An RNA-Seq case study added to User's Guide.

\item
getEAWP() now gets probe annotation from the expression rownames
  of an EList object, if no other probe annotation is available.

\item
topRomer() now ranks gene sets by secondary columns as well the
  primary criterion specified, to give a more meaningful ranking when
  the primary p-value is tied.

\item
change to p-value calculation in geneSetTest() when rank.only=FALSE.

\item
wilcoxGST() now accepts signed or unsigned test statistics.

\item
plotMA() now recognizes ElistRaw and EList objects appropriately.

\item
speed improvement for normalizeCyclicLoess when weights=NULL.

\item
default span for normalizeCyclicLoess increased from 0.4 to 0.7.

\item
bug fix for normalizeBetweenArrays() when object is an EListRaw and method="cyclicloess".  Previously this function was applying cyclicloess to the raw intensities, then logging.  Now it logs first, then applies cyclicloess.

\item
roast() and mroast() can now estimate and work with a trend on the prior variance.

\item
treat() now has an optional trend argument, same as eBayes().

\item
Weaver case study (Section 11.5) in User's Guide is updated and rewritten.

\item
new argument batch2 for removeBatchEffect().

\item
barcodeplot() and barcodeplot2() merged into one function.

\item
Data classes ElistRaw and Elist now described in the quick start section of the User's Guide.  Other minor updates to User's Guide.

\item
plotMDS is now an S3 generic function.  This allows MDS plots to be redrawn with new labels without needing to repeat the distance or scaling calculations.
New S4 class "MDS" to hold the multidimensional scaling information output from plotMDS.

\item
bug fix to avereps() for EList objects x when x$other is not empty.
  }
}
