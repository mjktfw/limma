\name{read.idat}
\alias{read.idat}

\title{Read Illumina expression data directly from IDAT files}

\description{Read Illumina BeadArray data from IDAT and manifest (.bgx) files for gene expression platforms.}

\usage{
read.idat(idatfiles, bgxfile, dateinfo=FALSE, tolerance=0)
}

\arguments{
  \item{idatfiles}{character vector specifying idat files to be read in.}
  \item{bgxfile}{character string specifying bead manifest file (.bgx) to be read in.}
  \item{dateinfo}{logical. Should date and software version info be read in?}
  \item{tolerance}{numeric. The number of probe ID discrepancies allowed between the manifest and a given idat file.} 
}

\details{
     Illumina's BeadScan/iScan software ouputs probe intensities in IDAT
     format (encrypted XML files) and probe info in a platform specific manifest file (.bgx).
     These files can be processed using the low-level functions \code{readIDAT} and \code{readBGX} 
     from the \code{illuminaio} package (Smith et al. 2013).  

     The \code{read.idat} function provides a convenient way to read these files.
     into R and store them in an \code{EListRaw-class} object (similar to \code{read.ilmn}, 
     which imports data output by Illumina's GenomeStudio software) that can be used 
     by downstream processing functions in \code{limma}.

     Probe types are indicated in the \code{Status} column of the \code{genes} 
     component of the \code{EListRaw-class} object.
}

\value{
  An \code{EListRaw-class} object with the following components:
  \item{E}{ numeric matrix of raw intensities.}
  \item{other}{ list containing matrices of \code{NumBeads} and \code{STDEV} for each probe.}
  \item{genes}{ data.frame of probe annotation.}
  \item{targets}{ data.frame of sample information.}
}

\references{
Smith ML, Baggerly KA, Bengtsson H, Ritchie ME, Hansen KD (2013). 
illuminaio: An open source IDAT parsing tool. \emph{F1000 Research} 2, 264.
\url{http://f1000research.com/articles/2-264/v1}
}

\author{Matt Ritchie}

\seealso{
     \code{read.ilmn} imports gene expression data output by GenomeStudio.

     \code{neqc} performs normexp by control background correction, log
     transformation and quantile between-array normalization for
     Illumina expression data.

     \code{propexpr} estimates the proportion of expressed probes in a microarray.
}

\examples{
\dontrun{
idatfiles = dir(pattern="idat")
bgxfile = dir(pattern="bgx")
data = read.idat(idatfiles, bgxfile)
propexpr(data)
datanorm = neqc(data)
}
}

\keyword{reading data}
\concept{illumina microarrays}
\concept{microarray data file}
