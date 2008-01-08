\name{multeq.rat}
\alias{multeq.rat}
\title{ Equivalence for ratios of means of multiple endpoints }
\description{
Performs equivalence tests and related confidence intervals for
ratios of two normal means of multiple endpoints.
}
\usage{
multeq.rat(data, grp, resp = NULL, margin.lo = NULL, margin.up = NULL, method = "single.step", FWER = 0.05)
}
\arguments{
  \item{data}{ a data frame containing response variables (endpoints) and the
               group variable as columns, the data must have exactly two treatment
               groups }
  \item{grp}{ the name of the group variable in " " }
  \item{resp}{ a vector of names of the response variables (endpoints) in " " }
  \item{margin.lo}{ a vector of relative lower margins under the null hypotheses
                    relating to the endpoints }
  \item{margin.up}{ a vector of relative upper margins under the null hypotheses
                    relating to the endpoints }
  \item{method}{ a character string: 
                   \item "step.up": method of Quan et al. (2001),
                   \item "single.step": Bonferroni-adjusted single-step 
                                        procedure}
  \item{FWER}{ a single numeric value specifying the familywise error rate to be
               controlled }
}
\details{
The object is to show equivalence for two treatment groups on multiple primary,
normally distributed response variables (endpoints). If margin.up is not given,
one-sided test are done for the alternative hypothesis that the ratios of the
means (first over second group) is larger than margin.lo. Analogously, same vice
versa. Only if both margin.lo and margin.up are given, a two-sided equivalence
tests for ratios is done. Bonferroni adjusted "two one-sided t-tests" (TOST)
and related simultaneous confidence intervals are used for method "single.step";
an extended version of the method of Quan et al. (2001) is applied for "step.up".
}
\value{
An object of class multeq.rat containing:
\item{estimate}{ a (named) vector of estimated ratios }
\item{test.stat}{ a (named) vector of the calculated test statistics }
\item{degr.fr}{ a single degree of freedom }
\item{p.value}{ a (named) vector of p-values adjusted for multiplicity }
\item{lower}{ a (named) vector of lower confidence limits }
\item{upper}{ a (named) vector of upper confidence limits }
}
\note{
Because related to the TOST method, the two-sided confidence intervals for
method="single.step" have simultaneous coverage probability (1-2alpha). The
intervals for method="step.up" are stepwise adjusted and only applicable for test
decisions, not for a simultaneous parameter estimation or comparing among each
other.
}
\references{
Quan et al. (2001): Assessmant of equivalence on multiple endpoints, Statistics in
Medicine 20, 3159-3173
}
\author{ Mario Hasler }
\seealso{ \code{\link{multeq.diff}} }
\examples{
data(clinic)

comp <- multeq.rat(data=clinic,grp="fact",method="step.up",margin.up=rep(1.25,5),
                   margin.lo=1/rep(1.25,5))
summary(comp)
}
\keyword{ htest }