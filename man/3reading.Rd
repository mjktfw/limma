\name{3.ReadingData}
\alias{3.Reading Data}
\title{Reading Microarray Data from Files}

\description{
This help page gives an overview of LIMMA functions used to read data into R from files.
}

\section{Reading Intensity Data}{
The first step in a microarray data analysis is to read into R the intensity data for each array provided by an image analysis program.
This is done using the function \code{\link{read.maimages}}.
This function is able to read output from SPOT, GenePix or QuantArray.

\code{\link{read.maimages}} optionally constructs quality weights for each spot using quality functions listed in \link{QualityWeights}.

\code{read.maimages} produces an \code{RGList} object and stores only the information required from each image analysis output file.
If you wish to read all the image analysis output files into R as individual data frames containing all the original columns, you may use \code{\link{read.series}}.
An \code{RGList} object can be extracted from the data frames at a later stage using the functions \code{\link{rg.spot}}, \code{\link{rg.genepix}} or \code{\link{rg.quantarray}}.

Another function, \code{\link{rg.series.spot}} is very similar to \code{\link{read.maimages}} with \code{source="spot"}.
This function will be removed in future versions of LIMMA.
}

\code{\link{read.maimages}} uses utility functions \code{\link{removeExt}} and \code{\link{read.matrix}}.

\section{Reading Gene Annotation}{

The function \code{\link{readGAL}} reads information from a GenePix Allocation List (gal) file.
It produces a data frame with known column names.

The function \code{\link{getLayout}} extracts from the gal-file data frame the print layout information for a spotted array.
The functions \code{\link{gridr}}, \code{\link{gridc}}, \code{\link{spotr}} and \code{\link{spotc}} use the extracted layout to compute grid positions and spot positions within each grid for each spot.

If each gene is printed more than once of the arrays, then \code{\link{uniquegenelist}} will remove duplicate names from the gal-file or gene list.

\code{\link{merge}} is designed to merge data sets with different gal-files provided that the gene library is the unchanged.
}

\author{Gordon Smyth}
\keyword{documentation}
