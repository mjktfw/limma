\name{limmanews}
\title{limma News}
\encoding{UTF-8}

\section{Version 3.22.0}{\itemize{
\item
New functions goana() and topGO() provide gene ontology analyses of differentially genes from a linear model fit.
The tests include the ability to adjust for gene length or abundance biases in differential expression detection, similar to the goseq package.

\item
Improvements to diffSplice.
diffSplice() now calculates Simes adjusted p-values for gene level inferences, in addition to the exon level t-tests and gene level F-tests.
topSplice() now has three ranking methods ("simes", "F" or "t"), with "simes" now becoming the default.
diffSplice() also has a new argument 'robust' giving access to robust empirical Bayes variance moderation.

\item
New function plotExons() to plot log-fold-changes by exon for a given gene.

\item
New function voomWithQualityWeights() allows users to estimate sample quality weights or allow for heteroscedasticity between treatment groups when doing an RNA-seq analysis.

\item
Improvement to arrayQualightyWeights().
It now has a new argument 'var.design' which allows users to model variability by treatment group or other covariates.

\item
Improved plotting for voomaByGroup().

\item
barcodeplot() can now plot different weights for different genes in the set.

\item
Improvements to roast() and mroast().
The directional (up and down) tests done by roast() now use both the original rotations and their opposite signs, effectively doubling the number of effective rotations for no additional computational cost.
The two-sided tests are now done explicitly by rotation instead of doubling the smallest one-sided p-value.
The two-sided p-value is now called "UpOrDown" in the roast() output.
Both functions now use a fast approximation to convert t-statistics into z-scores,
making the functions much faster when the number of rotations or the number of genes is large.
The contrast argument can now optionally be a character string giving a column name of the design matrix.

\item
zscoreT() can optionally use a fast approximation instead of the slower exact calculation.

\item
symbols2indices() renamed to ids2indices().

\item
Improvements to removeBatchEffect().
It can now take into account weights and other arguments that will affect the linear model fit.
It can now accept any arguments that would be acceptable for lmFit().
The behavior of removeBatchEffect() with design supplied has also changed so that it is now consistent with that of lmFit() when modelling batches as additive effects.
Previously batch adjustments were made only within the treatment levels defined by the design matrix.

\item
New function plotWithHighlights(), which is now used as the low-level function for plotMA() and plot() methods for limma data objects.

\item
The definition of the M and A axes for an MA-plot of single channel
  data is changed slightly.  Previously the A-axis was the average of
  all arrays in the dataset -- this has been definition since MA-plots
  were introduced for single channel data in April 2003.  Now an
  artificial array is formed by averaging all arrays other than the
  one to be plotted.  Then a mean-difference plot is formed from the
  specified array and the artificial array.  This change ensures the 
  specified and artificial arrays are computed from independent data,
  and ensures the MA-plot will reduce to a correct mean-difference
  plot when there are just two arrays in the dataset.

\item
plotMDS() can now optionally plot samples using symbols instead of
  text labels.
It no longer has a 'col' argument, which instead is handled by ....

\item
vennDiagram() now supports circles of different colors for any
  number of circles.  Previously this was supported only up to three
  sets.  

\item
getEAWP() will now find a weights matrix in an ExpressionSet object
  if it exists.

\item
update to helpMethods().

\item
Substantial updates to the two RNA-seq case studies in the User's Guide.
In both cases, the short read data has been realigned and resummarized.

\item
Improvements to many Rd files.
Many keyword entries have been revised.
Many usage and example lines been reformated to avoid over long lines.

\item
biocViews keywords updated.

\item
Subsetting columns of a MArrayLM object no longer subsets the design matrix.

\item
Bug fix for read.maimages: default value for 'quote' was not being
  set correctly for source="agilent.mean" or source="agilent.median".

\item
Bug fix to topTableF() and topTable().
The ordering of Amean values was sometimes incorrect when sorting by F-statistic and a lfc or p.value filter had been set.

\item
Bug fix to read.ilmn() when sep=",".
}}


\section{Version 3.20.0}{\itemize{

\item
New functions diffSplice(), topSplice() and plotSplice() provide functionality to analyse differential splicing using exon-level expression data from either microarrays or RNA-seq.

\item
New Pasilla case study added to User's Guide, demonstrating differential splicing analysis of RNA-seq data.

\item
new function weightedLowess() which fits a lowess curve with prior weights.
Unlike previous implementations of lowess or loess, the weights are used in calculating which neighbouring points to include in each local regression as well as in the local regression itself.

\item
weightedLoess() now becomes the default method used by loessFit() to fit the loess curve when there are weights.
The previous locfit and loess() methods are offered as options.

\item
linear model fit functions lm.series(), mrlm.series() and gls.series() no longer drop the dimensions of the components of the fitted object when there is just coefficient or just one gene.
Previously this was done inconsistently in some cases but not others.
Now the matrix components always keep dimensions.

\item
The functions lmFit(), eBayes() and tmixture.vector() now work even when there is just one gene (one row of data).

\item
New function subsetListOfArrays(), which is used to simplify the subsetting code for RGList, MAList, EList, EListRaw and MArrayLM objects.

\item
new function tricubeMovingAverage() for smoothing a time series.

\item
barcodeplot() has a new option to add enrichment worms to the plot, making use of tricubeMovingAverage().

\item
New plot() methods for RGList, MAList, EList and MArrayLM class objects.  In each case, this produces a similar result to plotMA().
When using plot() or plotMA() on an MArrayLM object, the column is now specified by the 'coef' argument instead of by 'array'.

\item
plotMA3by2() now works on single channel data objects as well as on MAList objects.

\item
New function read.idat() to read files from Illumina expression beadarrays in IDAT format.

\item
The ctrlpath argument of read.ilmn() now defaults to the same as path for regular probes.
This means that only one path setting is required if the regular and control probe profiles are in the same directory.

\item
read.ilmn() now sets the same probe IDs as rownames for both the expression matrix E and the annotation data.frame genes, providing that the probe IDs are unique.

\item
beadCountWeights() can now work with either probe-wise standard errors or probe-wise standard deviations.

\item
treat() has new arguments robust and winsor.tail.p which are passed through to robust empirical Bayes estimation.

\item
topTreat() now includes ... argument which is passed to topTable().

\item
topTable() with confint=TRUE now produces confidence intervals based on the t-distribution instead of on the normal distribution.
It also now accepts a numeric value for the confint argument to specify a confidence level other the default of 0.95.

\item
topTable() will now work on an MArrayLM fit object that is missing the lods component, for example as produced by treat().

\item
roast() and mroast() now permit array weights and observation weights to both be specified.

\item
camera(), roast() and mroast() now use getEAWP() to interpret the data object.  This means that they now work on any class of data object that lmFit() will.

\item
romer() now uses propTrueNull(method="lfdr") instead of convest().
This makes it substantially faster when the number of genes is large.

\item
genas() now uses fit$df.total from the MArrayLM object.
This prevents df.total from exceeding the total pooled residual df for the dataset.
The genas() results will change slightly for datasets for which df.prior was very lage.

\item
plotDensities() is now an S3 generic function with methods for RGList, MAList, EListRaw and EList objects.

\item
plotFB is now an S3 generic function with methods for RGList and EList data objects.

\item
New topic help pages 10GeneSetTests.Rd and 11RNAseq.Rd.
The page 10Other.Rd is deleted.
All topic help pages are now listed under 'See also' in the package introduction page accessible by ?limma.

\item
avereps() was never intended to be applied to RGList or EListRaw objects.
It now gives an error when applied to these objects instead of returning a matrix of questionable value.

\item
Bug fix: fitFDistRobustly() was failing when there were missing values or zero df values and covariate was NULL.

\item
Bug fix: vennDiagram() wasn't passing extra arguments (...) to plot() when the number of sets was greater than 3.

\item
Bug fix to topTreat().  Rownames were incorrectly ordered when p<1.

\item
bug fix to genas(), which was not handling vector df.prior correctly when the fit object was generated using robust=TRUE.

\item
bug fix to squeezeVar().
Previously there was an error when robust=TRUE and trend=FALSE and some of the estimated df.prior were infinite.

\item
bug fix to topTable() and topTableF() when sorting by F-statistic combined with p-value or lfc cutoffs.
}}


\section{Version 3.18.0}{\itemize{

\item
new function beadCountWeights() to estimate quantitative weights
from the bead counts for each probe for Illumina BeadArrays.

\item
New function contrastAsCoef(), which reforms a design matrix so
that one or more specified contrasts become coefficients.  This
function is called by roast().

\item
plotMA() is now an S3 generic function.

\item
The legend argument to plotMA() can now take a character value
giving the position to place the legend.  

\item
toptable(), topTable() and topTreat() now preserve the rownames of
the fit object, unless the fit object has duplicated names, in
which case the rownames are copied to the ID column.  Empty
rownames are replaced with 1:nrow(fit).

\item 
read.ilmn() no longer adds the probe IDs to the gene annotation
data.frame, leaving them instead as rownames of the expression
matrix.  It longer creates a targets file since the sample names
are already preserved as column names of the expression matrix.

\item
loessFit() now uses the locfit.raw in the locfit package when
weights are provided instead of loess in the stats package.  The
function now runs very efficiently even on very long data vectors.
The output results will change slightly when weights are provided.

\item
voom() now outputs lib.size as a column of targets instead of as a
separate component.

\item
cbind for EList and EListRaw objects now recognizes a design matrix
if it is present.

\item
plotMDS() now checks explicitly that there are at least 3 samples
to plot.

\item
normexp.fit.detection.p() now tolerates some non-monotonicity in
the detection p-pvalues as a function of expression.

\item
fitFDistRobustly() now uses a smoother for the smallest df.prior
values.  This may result in smaller tail values than before when
a group of input x values appear to be outliers but the largest
value is not individually a stand-out value.

\item
New merge methods for EList and EListRaw objects.

\item
topTable() and treat() now give more informative error messages
when the argument fit is not a valid MArrayLM object.

\item
roast() now calls mroast() if the index vector is a list.
Bug fix to mroast(), which had been ignoring the weights.

\item
Updates to genas() function.
argument chooseMethod renamed to subset and option "n" renamed to "all". 
Function now returns NA results and a message when no genes satisfy
the criterion for inclusion in the analysis.
Some editing of help page and streamlining of code.

\item
Roles of contributors now specified in author field of DESCRIPTION
file using standard codes.

\item
Additions and updates to references in the help pages.
Removed defunct Berkeley Press links to published Smyth (2004)
article in several Rd files.  Replacing with link to Preprint.
Added link to Phipson (2013) thesis in two Rd files.
Add Majewski et al reference to genas.Rd.
Add Phipson et al and Sartor et al references to squeezeVar.Rd.
Add Phipson et al reference to eBayes.Rd.
Update lmscFit and voom references.

\item
Update mammmary stem cell case study in User's Guide.  As well
as reflecting changes to read.ilmn() and topTable(), this now
demonstrates how to find signature genes for particular cell type.

\item
documentation about rownames and column names and the use of
rownames(fit) and colnames(fit) added to lmFit.Rd.

\item
improvements to help pages for data classes.

\item
Edits to normalizeBetweenArrays help page (i) to further clarify
which normalization methods are available for single-channel data
and which are available for two-color data and (ii) to give a cross
citation to the neqc() function for Illumina BeadChips.

\item
Edits to voomaByGroup help page.

\item
duplicateCorrelation() now uses the weights matrix when block is
set.  Previously the weights were ignored when block was used.

\item
Bug fix to subsetting for MArrayLM objects: the df.total component
was not being subsetted.

\item
bug fix to eBayes(robust=TRUE) when some of the df.prior values are
infinite.

\item
Bug fix to ebayes(), which was not passing the 'robust' argument
correctly on to squeezeVar().

\item
Bug fix to fitFDistRobustly(), which affected the estimated scale
when df2 is estimated to be Inf.
}}


\section{Version 3.16.0}{\itemize{

\item
New section in User's Guide on time course experiments with many time points.
The RNA-seq case study in User's Guide has also been revised.

\item
Improvements to various help pages including read.maimages.Rd, squeezeVar.Rd, fitFDist.Rd, trigammaInverse.Rd, normalizeRobustSpline.Rd, genas.Rd and roast.Rd.
Previously the meaning of source="agilent" was mis-stated in read.maimages.Rd.

\item
New robust method for estimating the empirical Bayes prior, called
by specifying robust=TRUE in the call to eBayes().  When this is
TRUE the output df.prior is now a vector instead of a scalar.

\item  
New function fitFDistRobustly() estimates the parameters of a
scaled F-distribution robustly using Winsorized values.  Outlier
observations receive smaller values for df.prior than non-outliers.
This permits robust methods for squeezeVar(), ebayes() and eBayes(),
all of which now have a new argument wins.tail.p to specify the
tail proportions for Winsorizing.

\item
fitFDist() now permits infinite values for the covariate.
It also gracefully handles cases where the covariate takes only a small number of distinct values.
Similarly for eBayes() and squeezeVar() that call fitFDist().  

\item
All the functions that perform gene set tests have been revised to make the input and output formats more consistent.

roast(), mroast() and camera() are now S3 generic functions, with methods for EList and MAList objects.

The order of arguments has been changed for roast(), mroast() and camera() so that the first argument is now y.

All functions that perform gene sets now use the argument 'index' to specify which genes are included in the test set.
Previously this argument was called 'iset' for roast() and romer() and 'indices' for camera().

camera() and mroast() now produce a data.frames.
Instead of separate up and down p-value columns, there is now a two-sided p-value and a column indicating direction of change.  
There are new columns giving FDR values and the number of genes in each set.
There is a new argument 'sort' to indicate whether output results should be sorted by p-value.

mroast() has a new argument 'weights' for observational weights, to bring it into line with roast(),

\item
vennDiagram() can now plot up to five sets (previously limited to three).

\item
genas() now optionally draws a plot in which ellipses are used to represent the technical and biological components of correlation.
It also now has the ability to automatically select which probes are used for the correlation analysis, and a new argument controls the method used for this selection.

\item
New options for the method argument of propTrueNull().

\item
New functions vooma() and voomaByGroup() for computing precision weights based on a mean-variance trend.
vooma() is similar to voom() but for microarray data instead of RNA-Seq.
voomaByGroup() allows different groups to have systematically different variances.

\item
New function predFCm() to compute predictive (shrunk) log fold changes.

\item
New function fitGammaIntercept() for estimating the intercept of a gamma glm with an offset.
Used by genas().

\item
New function zscoreHyper() for computing z-score equivalents of deviates from a hypergeometric distribution.

\item
New function qqf() for qq-plots relative to an F-distribution.

\item
normalizeWithinArrays() with method="robustspline" now longer requires the layout argument to be set.
The layout argument for normalizeRobustSpline() now defaults to a single print-tip group.

\item
fitFDist() now coerces degrees of freedom df1 below 1e-15 to zero.

\item
Due to changes in R, loessFit() no longer makes direct calls to foreign language code in the stats package, and instead calls R functions.
Unfortunately, this makes loessFit() about 25-30\% slower than previously when weights are used.

\item
Bug fix to read.maimages(), which was not accepting source="agilent.mean".

\item
Bug fix for contrasts.fit() when the covariance matrix of the coefficients (cov.coefficients) is not found in the fitted model object.
This situation doesn't arise using any of the standard limma analysis pipelines.

\item
Bug fix to lmscFit() when the residual df = 1.

\item
Bug fix to readTargets() to avoid warning message when targets$Label is used to set row names but targets$Label contains duplicated entries.

}}


\section{Version 3.14.0}{\itemize{

\item
limma license upgraded to GPL (>=2) instead of LGPL to match R itself.

\item
Many updates to the User's Guide. 
Sections have been added on reading single channel Agilent and Illumina data.
The chapter on experimental designs has been split into three chapters on single-channel, common reference and
two-color designs respectively.  The material on the fixed effect approach
  to technical replication has been deleted.
There are new sections on nested interactions for factorial designs
and on multi-level designs.

\item
The links to the Apoa1, Weaver and Bob1 datasets in the User's Guide have been updated
to help users download the data themselves if they wish to repeat the case study analyses.

\item
The help page for camera() now cites published paper Wu and Smyth (NAR, 2012).
In view of the results of this paper, the claim is no longer made on help page for geneSetTest()
  that genes might be treated as independent when the experimental
  units are genetically identical mice.

\item
Minor edits to CITATION file.

\item
New function propTrueNull() for fast estimation of the proportion of true null hypotheses from a vector of p-values.

\item
New function zscore() to compute z-score equivalents for deviates from any continuous distribution.
Includes the functionality of the older functions zscoreGamma() and zscoreT() as special cases.

\item
roast() now accepts observation level weights, through a new
  argument 'weights'.

\item
loessFit() now applies minimum and maximum bounds by default to avoid zero or
  infinite weights.  Equal weights are now treated as if the weights were NULL, even all
  zero weights, so that the lowess code is called instead of the loess code.

\item
When there are no weights, loessFit() now extracts residuals
directly from the C code output instead of computing in R.

\item
  fitFDist() now permits missing values for x or zero values for df1
  even when there is a covariate.  This means that squeezeVar() and
  eBayes() now work with trends even when not all the data values
  are informative.

\item
New argument 'file' for convest(), implementing edits contributed
  by Marcus Davy.  Arguments doplot and dereport renamed to 'plot'
  and 'report'.

\item
Two improvements for plotMDS().
It now coerces labels to be character, and now makes extra room on the plot when the text labels are wide.

\item
plotMDS() no longer gives an error when the requested number of top
  genes is greater than the total number of rows of data.

\item
Code speed-up for alias2SymbolTable()

\item
any(duplicated()) replaced by anyDuplicated() in several functions.

\item
Fix to voom() so that it computes weights correctly even when the
  design matrix is not of full rank.

\item
Bug fix for roast() when the fitted model has only one coefficient.

}}


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
