\name{normalizeBetweenArrays-generic}
\alias{normalizeBetweenArrays}
\title{Normalize Between Arrays - generic}

\description{
Normalizes expression intensities so that the intensities or log-ratios have similar distributions across a series of arrays.
}

\usage{
normalizeBetweenArrays(object,method)
}

\arguments{
  \item{object}{an object, usually an \code{\link[limma:malist]{MAList}}, containing expression ratios for a series of arrays}
  \item{method}{character string specifying the normalization method to be used.
  Choices are \code{"none"}, \code{"scale"} or \code{"quantile"}.
  A partial string sufficient to uniquely identify the choice is permitted.
  If no method is recognized then unnormalized values are returned.}
}

\details{
\code{normalizeWithinArrays} normalizes expression values to make intensities consistent within each array.
\code{normalizeBetweenArrays} normalizes expression values to achieve consistency between arrays.

The scale normalization method was proposed by Yang et al (2001, 2002) and is further explained by Smyth and Speed (2003).
The idea is simply to scale the log-ratios to have the same median-abolute-deviation (MAD) across arrays.
This idea has also been implemented by the \code{maNormScale} function in the marrayNorm package.
The implementation here is slightly different in that the MAD scale estimator is replaced with the median-absolute-value and the A-values are normalized as well as the M-values.

Quantile normalization was proposed by Bolstad et al (2003) for Affymetrix-style single-channel arrays and by Yang and Thorne (2003) for two-color cDNA arrays.
}

\section{Methods}{
\item{\code{object=MAList, method=="scale"}}{
M-values are scaled to have the same median-absolute values across arrays, A-values are scaled to have the same medians across arrays
}
\item{\code{object=MAList, method=="quantile"}}{
M- and A-values are normalized so that the single channel intensities R and G have the same empirical distribution across channels and across arrays
}
\item{\code{object=list}}{
\code{object} is converted into an \code{MAList} if not one already
}
\item{\code{object=matrix, method=="scale"}}{
If \code{object} is a plain matrix, then its columns are normalized to have equal median-absolute-values.
}
\item{\code{object=matrix, method=="quantile"}}{
If \code{object} is a plain matrix, then its columns are normalized to have equal empirical distributions.
}
}

\value{
If \code{object} is a matrix then \code{normalizeBetweenArrays} produces a matrix of the same size.
Otherwise, \code{normalizeBetweenArrays} produces an \code{\link[limma:malist]{MAList}} object.
}

\author{Gordon Smyth}

\references{
Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003), A comparison of normalization methods for high density oligonucleotide array data based on bias and variance. \emph{Bioinformatics} \bold{19}, 185-193.

Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray data. In: \emph{METHODS: Selecting Candidate Genes from DNA Array Screens: Application to Neuroscience}, D. Carter (ed.). To appear.

Yang, Y. H., Dudoit, S., Luu, P., and Speed, T. P. (2001). Normalization for cDNA microarray data. In \emph{Microarrays: Optical Technologies and Informatics}, M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), Proceedings of SPIE, Volume 4266, pp. 141-152. 

Yang, Y. H., Dudoit, S., Luu, P., Lin, D. M., Peng, V., Ngai, J., and Speed, T. P. (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research} \bold{30}(4):e15.

Yang, Y. H., and Thorne, N. P. (2003). Normalization for two-color cDNA microarray data.
In: D. R. Goldstein (ed.), \emph{Science and Statistics: A Festschrift for Terry Speed}, IMS Lecture Notes - Monograph Series, Volume 40, pp. 403-418.
}

\seealso{
  An overview of LIMMA functions for normalization is given in \link{4.Normalization}.

  See also \code{\link[marrayNorm]{maNormScale}} in the marrayNorm package and \code{\link[affy:normalize-methods]{normalize}} in the affy package.
}

\examples{
library(sma)
data(MouseArray)
MA <- normalizeWithinArrays(mouse.data, mouse.setup)
plot.scale.box(MA$M)
MA <- normalizeBetweenArrays(MA,method="scale")
plot.scale.box(MA$M)

#  One can get the same results using the matrix method:
M <- normalizeBetweenArrays(MA$M,method="scale")
plot.scale.box(M)
}

\keyword{models}
\keyword{multivariate}
