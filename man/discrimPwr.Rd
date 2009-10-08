\name{discrimPwr}
\alias{discrimPwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensory discrimination power analysis}
\description{
  Computes the power of the hypothesis test of no sensory
  difference for any one of four methods: 2-AFC, 3-AFC, duotrio and
  triangle tests given the underlying sensory difference delta, the type I
  test level and the sample size.
}
\usage{
discrimPwr(delta, sample.size, alpha = 0.05,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
           pd0 = 0, type = c("difference", "similarity"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{the underlying sensory difference (non-negative)}
  \item{sample.size}{the sample size (a positive integer)}
  \item{alpha}{the type I level of the test (must be between zero and one)}
  \item{method}{the discrimination test protocol. Four allowed values:
    \code{"twoAFC"}, \code{"threeAFC"}, \code{"duotrio"},
    \code{"triangle"}}
  \item{pd0}{the proportion of discriminators in the population of interest}
  \item{type}{the type of test}
}
\details{
  The power of the standard one-tailed difference test of "no
  difference" is obtained with \code{pd0 = 0}.

  The probability under the null hypothesis is
  given by \code{pd0 + p0 * (1 - pd0)} where \code{p0} is the guessing
  probability defined by the \code{method} argument.

  The function uses one of the dedicated binomial families.
}
\value{
  The power.
}
\references{Brockhoff, P.B. and Christensen, R.H.B (2010). Thurstonian
  models for sensory discrimination tests as generalized linear models.
  Food Quality and Preference, 21, pp. 330-338.}
\author{Rune Haubo B Christensen and Per Bruun Brockhoff}

\seealso{\code{\link{triangle}}, \code{\link{twoAFC}},
  \code{\link{threeAFC}}, \code{\link{duotrio}},
  \code{\link{discrim}}, \code{\link{discrimSim}},
  \code{\link{AnotA}}, \code{\link{discrimSS}},
  \code{\link{samediff}}, \code{\link{findcr}}
}

\examples{
## Finding the power of a discrimination test with a sensory delta of 1,
## a sample of size 30 and a type I level of .05:
discrimPwr(1, 30, 0.05, "twoAFC")
discrimPwr(1, 30, 0.05, "threeAFC")
discrimPwr(1, 30, 0.05, "duotrio")
discrimPwr(1, 30, 0.05, "triangle")

## A similarity example:
discrimPwr(.2, 100, method = "triangle", pd0 = .2, type = "simil")

}
\keyword{models}% at least one, from doc/KEYWORDS
