\name{3.ReadingData}
\alias{3.ReadingData}
\title{Reading Microarray Data from Files}

\description{
This help page gives an overview of LIMMA functions used to read data into R from files.
}

\section{Reading Intensity Data}{
The first step in a microarray data analysis is to read into R the intensity data for each array provided by an image analysis program.
This is done using the function \code{\link{read.maimages}}.

\code{\link{read.maimages}} optionally constructs quality weights for each spot using quality functions listed in \link{QualityWeights}.

\code{read.maimages} produces an \code{RGList} object and stores only the information required from each image analysis output file.
If you wish to read all the image analysis output files into R as individual data frames containing all the original columns, you may use \code{\link{read.series}}.
An \code{RGList} object can be extracted from the data frames at a later stage using the functions \code{\link{rg.spot}}, \code{\link{rg.genepix}} or \code{\link{rg.quantarray}}.

Another function, \code{\link{rg.series.spot}} is very similar to \code{\link{read.maimages}} with \code{source="spot"}.
This function will be removed in future versions of LIMMA.
}

\code{\link{read.maimages}} uses utility functions \code{\link{removeExt}}, \code{\link{read.matrix}}, \code{\link{read.imagene}} and \code{\link{getImageneHeaders}}.

\section{Reading Gene Annotation}{

Many image analysis program provide gene IDs as columns in the image analysis output files, for example ArrayVision, Imagene and the Stanford Microarray Database.
In other cases you may have gene name and annoation information in a separate file.
The function \code{\link{readGAL}} reads information from a GenePix Allocation List (gal) file.
It produces a data frame with known column names.
If the gene names consist of a short name followed by annotation information, then \code{\link{splitName}} may be used to separate the name and annotation information into separate vectors.

The function \code{\link{getLayout}} extracts from the gal-file data frame the print layout information for a spotted array.
The functions \code{\link{gridr}}, \code{\link{gridc}}, \code{\link{spotr}} and \code{\link{spotc}} use the extracted layout to compute grid positions and spot positions within each grid for each spot.
The function \code{\link{printorder}} calculates the printorder, plate number and plate row and column position for each spot given information about the printing process.

If each gene is printed more than once of the arrays, then \code{\link{uniquegenelist}} will remove duplicate names from the gal-file or gene list.

\code{\link[limma:cbind]{cbind}} allows different \code{RGList} or \code{MAList} objects to be combined assuming the layout of the arrays to be the same.
\code{\link[limma:merge]{merge}} can combine data even when the order of the genes on the arrays has changed.
\code{merge} uses utility function \code{\link{makeUnique}}.
}

\author{Gordon Smyth}
\keyword{documentation}
