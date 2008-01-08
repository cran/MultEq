\name{print.multeq.diff}
\alias{print.multeq.diff}
\title{ Print out of the results of multeq.diff }
\description{
A short print out of the results of multeq.diff.
}
\usage{
\method{print}{multeq.diff}(x, digits = 4, ...)
}
\arguments{
  \item{x}{ an object of class "multeq.diff" as obtained by calling multeq.diff }
  \item{digits}{ digits for rounding the results }
  \item{\dots}{ arguments to be passed to print }
}
\value{
A print out containing the margins, estimates, confidence intervals, and 
p.values computed by multeq.diff.
}
\author{ Mario Hasler }
\seealso{ \code{\link{print.multeq.rat}} }
\keyword{ print }
