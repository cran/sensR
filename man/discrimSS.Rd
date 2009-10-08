\name{discrimSS}
\alias{discrimSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sensory discrimination sample size calculation}
\description{
   Computes the sample size for any one of four methods: 2-AFC, 3-AFC,
   duotrio and triangle tests given
   the underlying sensory difference d, the type 1 test level and the
   required power.
}
\usage{
discrimSS(delta, power, alpha = 0.05,
          method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
          pd0 = 0, type = c("difference", "similarity"), start = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{the underlying sensory difference (larger than zero)}
  \item{power}{the wanted power (between zero and one)}
  \item{alpha}{the type 1 level of the test (between zero and one)}
  \item{method}{the discrimination protocol. Four allowed values:
    "twoAFC", "threeAFC", "duotrio", "triangle"}
  \item{pd0}{the proportion of discriminators in the population of interest}
  \item{type}{the type of test to be conducted}
  \item{start}{lower bound on the sample size. Specifying \code{start} at
  a value close to, but lower than the sample size will save
  computational time. If a too high value is given, this value is
  returned and the user can try a lower value for \code{start}.}
}
\details{
  The sample size of the standard one-tailed difference test of "no
  difference" is obtained with \code{pd0 = 0}.

  The probability under the null hypothesis is
  given by \code{pd0 + p0 * (1 - pd0)} where \code{p0} is the guessing
  probability defined by the \code{method} argument.

  The function uses one of the dedicated binomial families and the
  \code{discrimPwr} function
}
\value{
  The sample size
}
\references{Brockhoff, P.B. and Christensen, R.H.B (2010). Thurstonian
  models for sensory discrimination tests as generalized linear models.
  Food Quality and Preference, 21, pp. 330-338.}
\author{Per Bruun Brockhoff and Rune Haubo B Christensen}

\seealso{ \code{\link{triangle}}, \code{\link{twoAFC}},
  \code{\link{threeAFC}}, \code{\link{duotrio}},
  \code{\link{discrim}}, \code{\link{discrimSim}},
  \code{\link{AnotA}}, \code{\link{discrimPwr}},
  \code{\link{samediff}}, \code{\link{findcr}}}

\examples{
## Finding the necessary sample size:
discrimSS(1, 0.9, 0.05, "twoAFC")
discrimSS(1, 0.9, 0.05, "threeAFC")
discrimSS(1, 0.9, 0.05, "duotrio")
discrimSS(1, 0.9, 0.05, "triangle")

## A similarity example:
discrimSS(.2, .80, method = "triangle", pd0 = .2, type = "simil")


}
\keyword{models}% at least one, from doc/KEYWORDS
