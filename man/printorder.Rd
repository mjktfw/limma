\name{printorder}
\alias{printorder}
\title{Identify Order in which Spots were Printed}
\description{
Identify order in which spots were printed and the 384-well plate from which they were printed.
}
\usage{
printorder(layout,nwells=384,ntips=layout$ngrid.r*layout$ngrid.c)
}
\arguments{
  \item{layout}{list with the components \code{ngrid.r}, \code{ngrid.c}, \code{nspot.r} and \code{nspot.c}}
  \item{nwells}{number of wells in plate containing DNA to be printed}
  \item{ntips}{actual number of tips of the printer-head}
}
\details{
In most cases the printer-head contains the \code{layout$ngrid.r} times \code{layout$ngrid.c} tips and the array is printed using \code{layout$nspot.r} times \code{layout$npot.c} dips of the head.

In some cases a smaller number of physical tips is used and the total number of grids is built up by effectively printing two or more sub-arrays on the same slide.
In this case the number of grids should be a multiple of the number of tips.
}
\value{
List with components
\item{printorder}{numeric vector giving printorder of each spot}
\item{plate}{numeric vector giving plate number of each spot}
}
\author{Gordon Smyth}
\keyword{IO}
