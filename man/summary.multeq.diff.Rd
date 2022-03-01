\name{summary.multeq.diff}
\alias{summary.multeq.diff}
\title{ Summary function for multeq.diff }
\description{
A detailed print out of the results of multeq.diff.
}
\usage{
\method{summary}{multeq.diff}(object, digits = 4, ...)
}
\arguments{
  \item{object}{ an object of class "multeq.diff" as obtained by calling multeq.diff }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the margins, degrees of freedom, estimates, test
statistics, confidence intervals, and p.values computed by multeq.diff.
}
\author{ Mario Hasler }
\seealso{ \code{\link{summary.multeq.rat}} }
\keyword{ print }