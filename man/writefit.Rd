\name{write.fit}
\alias{write.fit}

\title{Write MArrayLM Pbject to a File}

\description{
Write a microarray linear model fit a file.
}

\usage{
write.fit(fit, clas=NULL, file, digits=2)
}

\arguments{
\item{fit}{object of class \code{MArrayLM} containing the results of a linear model fit}
\item{clas}{object of class \code{classification}}
\item{file}{character string giving name of file}
\item{digits}{integer indicating precision to be used}
}

\value{
No value is produced but a file is written.
}

\details{
This function writes a tab-delimited text file containing for each gene (1) the average log-intensity, (2) the log-ratios, (3) moderated t-statistics, (4) F-statistic if available, (5) classification if available and (6) gene names and annotation.
}

\author{Gordon Smyth}
\seealso{
\code{\link[base]{write}} in the base library.

An overview of linear model functions in limma is given by \link{5.LinearModels}.
}
\keyword{IO}
