\name{limmanews}
\title{limma News}
\encoding{UTF-8}

\section{Version 3.12.0}{\itemize{

\item
read.maimages() with source="agilent" now reads median foreground
estimates instead of mean foreground.  New option source=
"agilent.mean" preserves earlier meaning of source="agilent".

\item
Agilent single-channel case study added to User's Guide.

\item
removeBatchEffect() now corrects for continuous covariates
as well as qualitative factors.

\item
new function camera() performs competitive
gene set tests while adjusting for inter-gene correlation.

\item
new function interGeneCorrelation() estimates the average
intergene correlation for a set of genes.

\item
columns in output from roast() have been re-ordered.

\item
arguments 'selected' and 'selected2' renamed to 'index'
and 'index2' in functions barcodeplot(), geneSetTest()
and wilcoxGST().

\item
default labels for barcodeplot() are now somewhat more
explicit.

\item
new function rankSumTestWithCorrelation extends the
Wilcoxon-Mann-Whitney test to allow for correlation between
cases in one of the groups.  geneSetTest() now calls this
function instead of wilcox.test, with a consequence improvement
in speed.

\item
The lfc (log-fold-change) cutoff argument of topTable() is now
applied to the minimum absolute logFC when ranking by F-statistic.
Previously lfc was only used when ranking by t-statistic.

\item
new methods "fast" and "affy" for normalizeCyclicLoess(),
with "fast" becoming the default method.
New argument 'cyclic.method' for normalizeBetweenArrays() gives
access to the different cyclic loess methods.

\item
There were problems with using the argument gene.weights in
mroast().  This argument is now permitted to be of the same
length as the number of probes in the data set.  It is then
automatically subsetted for each gene set.

\item
mroast() now uses mid-p-values by default when adjusting for multiple testing.

\item
neqc(), nec() and normexp.fit.control() now give user-friendly error messages when no negative control probes or no regular probes are found.

}}

\section{Version 3.10.0}{\itemize{

\item
New function voom() allows RNA-Seq experiments to be analysed using the standard limma pipeline.
An RNA-Seq case study is added to User's Guide.

\item
treat(), roast() and mroast() can now estimate and work with a trend on the prior variance, bringing them into line with eBayes().

\item
barcodeplot() and barcodeplot2() merged into one function.

\item
removeBatchEffect() can now correct for two batch factors.

\item
plotMDS is now an S3 generic function.  This allows MDS plots to be redrawn with new labels without needing to repeat the distance or scaling calculations.
New S4 class "MDS" to hold the multidimensional scaling information output from plotMDS.

\item
getEAWP() now gets probe annotation from the expression rownames of an EList object, if no other probe annotation is available.

\item
topRomer() now ranks gene sets by secondary columns as well the primary criterion specified, to give a more meaningful ranking when the p-values are tied.

\item
wilcoxGST() now accepts signed or unsigned test statistics.
Change to p-value calculation in geneSetTest() when rank.only=FALSE to avoid zero p-values and follow recommendation of Phipson and Smyth (SAGMB, 2010).

\item
plotMA() now recognizes ElistRaw and EList objects appropriately.

\item
Default span for normalizeCyclicLoess increased from 0.4 to 0.7.
Speed improved when weights=NULL.

\item
Weaver case study (Section 11.5) in User's Guide is updated and rewritten.
Data classes ElistRaw and Elist now described in the quick start section of the User's Guide.  Other minor updates to User's Guide.

\item
Bug fix for normalizeBetweenArrays() when object is an EListRaw and method="cyclicloess".  Previously this function was applying cyclicloess to the raw intensities, then logging.  Now it logs first, then applies cyclicloess.

\item
Bug fix to avereps() for EList objects x when x$other is not empty.

}}
