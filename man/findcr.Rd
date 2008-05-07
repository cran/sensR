\name{findcr}
\alias{findcr}
\title{Find the critical value of a one-tailed binomial test}
\description{Finds the critical value in a one-tailed binomial test}
\usage{
findcr(sample.size, alpha = .05, p0 = .5)
}
\arguments{
  \item{sample.size}{the sample size of the binomial test (must be a positve integer)}
  \item{alpha}{the type I error-level of the test (must be between zero and one)}
  \item{p0}{the binomial probability under the null-hypothesis (must be
    between zero and one); 1/2 for the duotrio and twoAFC tests and 1/3
    for the triangle and threeAFC tests}
}
\value{
The critical value in a one-tailed binomial test, that is, the smallest integer such that 
the null hypothesis binomial probability of being larger than or equal to this number is 
smaller than or equal to the type I error-level of the test
}
\author{Rune Haubo B Christensen and Per Bruun Brockhoff}

\seealso{\code{\link{triangle}}, \code{\link{twoAFC}},
  \code{\link{threeAFC}}, \code{\link{duotrio}},
  \code{\link{discrim}}, \code{\link{discrimPwr}},
  \code{\link{discrimSim}}, \code{\link{AnotA}}
  \code{\link{discrimSS}}, \code{\link{samediff}}}

\examples{
## Find the critical value for a triangle test for the level 0.05 test 
## with 25 subjects:
findcr(sample.size = 25, , p0 = 1/3) 
}
\keyword{models}% at least one, from doc/KEYWORDS

