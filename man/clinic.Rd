\name{clinic}
\alias{clinic}
\docType{data}
\title{ Body measurements in a clinical study }
\description{
Measurements on six parts of patients' bodies in a clinical study for two competing
treatments.
}
\usage{data(clinic)}
\format{
  A data frame with 30 observations on the following 6 variables.
  \describe{
    \item{\code{fact}}{a factor with levels \code{1} \code{2}, specifying the treatment
      groups }
    \item{\code{var1}}{numeric vectors containing measurements on a first part of
      patients' bodies}
    \item{\code{var2}}{numeric vectors containing measurements on a second part of
      patients' bodies}
    \item{\code{var3}}{numeric vectors containing measurements on a third part of
      patients' bodies}
    \item{\code{var4}}{numeric vectors containing measurements on a fourth part of
      patients' bodies}
    \item{\code{var5}}{numeric vectors containing measurements on a fifth part of
      patients' bodies}
  }
}
\source{
L"auter, and Kropf, (1998): Exact stable multivariate tests for application in clinical
research. Joint statistical meeting Dallas (USA), conference proceeedings, group 1
}
\examples{
library(MultEq)

data(clinic)
plot(clinic[,-1])
}
\keyword{datasets}
