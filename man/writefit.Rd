\name{write.fit}
\alias{write.fit}

\title{Write MArrayLM Object to a File}

\description{
Write a microarray linear model fit to a file.
}

\usage{
write.fit(fit, results=NULL, file, digits=3, adjust="none", method="separate",
          F.adjust="none", sep="\t", \dots)
}

\arguments{
\item{fit}{object of class \code{MArrayLM} containing the results of a linear model fit}
\item{results}{object of class \code{TestResults}}
\item{file}{character string giving name of file}
\item{digits}{integer indicating precision to be used}
\item{adjust}{character string specifying multiple-testing adjustment method for the t-statistic P-values, e.g., \code{"BH"}. See \code{\link[stats]{p.adjust}} for the available options. If \code{NULL} or \code{"none"} then the P-values are not adjusted.}
\item{method}{character string, should the P-value adjustment be \code{"global"} or \code{"separate"} for each contrast.}
\item{F.adjust}{character string specifying adjustment method for the F-statistic P-values.}
\item{sep}{the field separator string. Values in the output file will be separated by this string.}
\item{\dots}{other arguments are passed to \code{write.table}}
}

\value{
No value is produced but a file is written to the current working directory.
}

\details{
This function writes a tab-delimited text file containing for each gene (1) the average log-intensity, (2) the log-ratios, (3) moderated t-statistics, (4) t-statistic P-values, (5) F-statistic if available, (6) F-statistic P-values if available, (7) classification if available and (8) gene names and annotation.
}

\author{Gordon Smyth}
\seealso{
\code{\link[base]{write}} in the base library.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\keyword{IO}
