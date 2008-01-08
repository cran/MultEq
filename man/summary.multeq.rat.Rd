\name{summary.multeq.rat}
\alias{summary.multeq.rat}
\title{ Summary function for multeq.rat }
\description{
A detailed print out of the results of multeq.rat.
}
\usage{
\method{summary}{multeq.rat}(object, digits = 4, ...)
}
\arguments{
  \item{object}{ an object of class "multeq.rat" as obtained by calling multeq.rat }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the margins, degrees of freedom, estimates, test
statistics, confidence intervals, and p.values computed by multeq.rat.
}
\author{ Mario Hasler }
\seealso{ \code{\link{summary.multeq.diff}} }
\keyword{ print }