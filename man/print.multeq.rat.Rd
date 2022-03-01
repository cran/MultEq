\name{print.multeq.rat}
\alias{print.multeq.rat}
\title{ Print out of the results of multeq.rat }
\description{
A short print out of the results of multeq.rat.
}
\usage{
\method{print}{multeq.rat}(x, digits = 4, ...)
}
\arguments{
  \item{x}{ an object of class "multeq.rat" as obtained by calling multeq.rat }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the margins, estimates, confidence intervals, and 
p.values computed by multeq.rat.
}
\author{ Mario Hasler }
\seealso{ \code{\link{print.multeq.diff}} }
\keyword{ print }
