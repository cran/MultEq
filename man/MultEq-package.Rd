\name{MultEq-package}
\alias{MultEq-package}
\alias{MultEq}
\docType{package}
\title{
Equivalence for multiple endpoints
}
\description{
The package provides tests and confidence intervals for comparing two
treatments when there is more than one primary response variable
(endpoint). The step-up procedure of Quan et al. (2001) is both
applied for differences and extended to ratios of means of normally
distributed data with equal group variances. A related single-step
procedure is also available.
}
\details{
\tabular{ll}{
Package: \tab MultEq\cr
Type: \tab Package\cr
Version: \tab 2.3\cr
Date: \tab 2011-10-20\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
\itemize{
  \item{multeq.diff}{Equivalence tests and related confidence intervals for
    differences of normal means of multiple endpoints}
  \item{multeq.rat}{Equivalence tests and related confidence intervals for
    ratios of normal means of multiple endpoints}
  \item{clinic}{Data set of body measurements in a clinical study} }
}
\author{
Mario Hasler

Maintainer: Mario Hasler <hasler@email.uni-kiel.de>
}
\references{
Quan et al. (2001): Assessmant of equivalence on multiple endpoints,
Statistics in Medicine 20, 3159-3173
}
\keyword{ package }
\examples{
data(clinic)

comp <- multeq.diff(data=clinic,grp="fact",method="step.up",margin.up=rep(0.6,5),
                    margin.lo=-rep(0.6,5))
summary(comp)
}
