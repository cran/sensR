\name{discrimSS}
\alias{discrimSS}
\alias{d.primeSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensory discrimination sample size calculation}
\description{
  Computes the sample size for a difference or similarity test for a
  sensory discrimination experiment using the binomial distribution. The
  function returns the smallest sample size for the test that gives the
  at least the desired target-power for the one-sided exact binomial
  test. A normal approximation to the sample size / power may used.
  \code{d.primeSS} is a convenience function that calls
  \code{discrimSS} but has arguments in terms of d-prime rather than
  pd, the probability of discrimination.
}
\usage{
discrimSS(pdA, pd0 = 0, target.power = 0.90, alpha = 0.05,
          pGuess = 1/2, test = c("difference", "similarity"),
          statistic = c("exact", "normal"))

d.primeSS(d.primeA, d.prime0 = 0, target.power = 0.90, alpha = 0.05,
          method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
          test = c("difference", "similarity"),
          statistic = c("exact", "normal")) 
}
\arguments{
  \item{pdA}{the probability of discrimination for the model under the
    alternative hypothesis; scalar between zero and one}
  \item{d.primeA}{d-prime for the model under the
    alternative hypothesis; non-negative numerical scalar}
  \item{pd0}{the probability of discrimination under the
    null hypothesis; scalar between zero and one}
  \item{d.prime0}{d-prime under the null hypothesis; non-negative
    numerical scalar}
  \item{target.power}{the desired power for the test}
  \item{alpha}{the type I level of the test; scalar between zero and
    one}
  \item{method}{the discrimination protocol for which the sample size
    should be computed} 
  \item{pGuess}{the guessing probability for the discrimination
    protocol, e.g. 1/2 for duo-trio and 2-AFC, and 1/3 for
    triangle and 3-AFC; scalar between zero and one}
  \item{test}{the type of one-sided binomial test (direction of the
    alternative hypothesis): "difference" corresponds "greater" and
    "similarity" corresponds to "less"}
  \item{statistic}{should sample size calculation be based on the 'exact'
    binomial test or the normal approximation to this?}
}
\details{
  The sample size of the standard one-tailed difference test where the
  null hypothesis of "no difference" is obtained with \code{pd0 = 0}.

  The probability under the null hypothesis is
  given by \code{pd0 + pg * (1 - pd0)} where \code{pg} is the guessing
  probability \code{pGuess}. Similarly, the probability of the
  alternative hypothesis is given by \code{pdA + pg * (1 - pdA)}

  The sample size for the 'exact' test, the function uses 
  the normal approximation to get a starting value for an iterative
  search. If the sample size based on the normal approximation is larger
  than 10,000, the function returns the normal approximation and issues
  a warning.
}
\value{
  The sample size; a scalar integer.
}
\references{Brockhoff, P.B. and Christensen, R.H.B (2010). Thurstonian
  models for sensory discrimination tests as generalized linear models.
  Food Quality and Preference, 21, pp. 330-338.}
\author{Per Bruun Brockhoff and Rune Haubo B Christensen}

\seealso{\code{\link{AnotA}}, \code{\link{discrimPwr}},
  \code{\link{samediff}}, \code{\link{findcr}}} 

\examples{
## Finding the smallest necessary sample size:
discrimSS(pdA = 0.5, pd0 = 0, target.power = 0.80, alpha = 0.05,
   pGuess = 1/2, test = "difference", statistic = "exact")

## Give identical results:
pd <- coef(rescale(d.prime = 1, method = "twoAFC"))$pd
discrimSS(pdA = pd, pd0 = 0, target.power = 0.90, alpha = 0.05,
   pGuess = 1/2, test = "difference", statistic = "exact")
d.primeSS(1, target.power = 0.90, method = "twoAFC")

## A similarity example:
discrimSS(pdA = 0.1, pd0 = 0.2, target.power = 0.80, alpha = 0.05,
   pGuess = 1/2, test = "similarity", statistic = "exact")

}
\keyword{models}% at least one, from doc/KEYWORDS
