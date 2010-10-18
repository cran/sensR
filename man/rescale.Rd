\name{rescale}
\alias{rescale}
%% \alias{print.rescale}
\alias{psyfun}
\alias{psyinv}
\alias{psyderiv}
\alias{pc2pd}
\alias{pd2pc}
\title{Transform or rescale between Pc, Pd and d-prime for sensory
  discrimination protocols}
\description{Transforms or rescales estimates and optionally standard
  errors between the three levels at which a sensory difference is
  measured: Pc (proportion of correct answers), Pd (proportion of
  discriminators) and d-prime. \code{rescale} is the main function and
  only one of Pc, Pd or d-prime should be given as argument --- values
  for the remaining two scales will be computed.

  A number of auxiliary functions are also provided:
  \describe{
    \item{\code{psyfun}}{implements the psychometric functions and maps
      from d-prime to Pc}
    \item{\code{psyinv}}{implements the inverse psychometric functions
      and maps from Pc to d-prime}
    \item{\code{psyderiv}}{implements the derivative of the psychometric
      functions}
    \item{\code{pc2pd}}{maps from Pc to Pd}
    \item{\code{pd2pc}}{maps from Pd to Pc}
    }
} 
\usage{

rescale(Pc, Pd, d.prime, std.err, 
        method = c("duotrio", "threeAFC", "twoAFC", "triangle"))

psyfun(d.prime, method = c("duotrio", "threeAFC", "twoAFC", "triangle"))

psyinv(Pc, method = c("duotrio", "threeAFC", "twoAFC", "triangle"))

psyderiv(d.prime, method = c("duotrio", "threeAFC", "twoAFC", "triangle"))

pc2pd(Pc, Pguess)

pd2pc(Pd, Pguess)

}
\arguments{
  \item{Pc}{the proportion of correct answers; a numerical vector
    between 0 and 1}
  \item{Pd}{the proportion of discriminators; a numerical vector between
    0 and 1} 
  \item{d.prime}{the sensory difference on the d-prime scale; a
    non-negative numerical vector.} 
  \item{std.err}{optional numerical vector of standard errors of the
    same length as the either of \code{Pc}, \code{Pd} or
    \code{d.prime}. Negative values are not allowed, but values may be
    \code{NA}}  
  \item{method}{the sensory discrimination protocol for which the
    results should apply}
  \item{Pguess}{the guessing probability implied by the protocol; a
    numeric scalar between 0 and 1}
}
\details{
  The \code{rescale} function is based on the fact that once the
  protocol and one of 
  Pc, Pd and d-prime is known, the other two can be computed. The same
  applies to the standard errors of these parameters.

  Standard errors are optional, but if they are supplied, the length of
  the \code{std.err} argument has to match the length of \code{Pc},
  \code{Pd} or \code{d.prime} whichever is given.

  A \code{print} method is implemented for \code{rescale} objects.
}
\value{
  For \code{rescale} an object of class \code{rescale} with elements
  \item{coefficients}{a \code{data.frame} with values of \code{Pc},
    \code{Pd} and \code{d.prime} corresponding to the input}
  \item{std.err}{if standard errors are given trough the \code{std.err}
    argument a \code{data.frame} of the same size and shape as
    \code{coefficients} with standard errors. Otherwise missing.}
  \item{method}{the sensory discrimination protocol for which the
    results apply}

  For \code{psyfun}, \code{psyinv}, \code{psyderiv}, \code{pc2pd} and
  \code{pd2pc} a numerical vector of the same length as the first
  argument with appropriate contents.
}

\author{Rune Haubo B Christensen}

\examples{

## suppose 15 out of 20 are observed in a duo-trio experiment, then
## the estimated probability of correct a answer is
(Pc <- 15/20)
## The standard error of this estimate is
(se.Pc <- sqrt(Pc * (1 - Pc) / 20))
## The corresponding estimate of proportion of discriminators (Pd) and
## d-prime with associated standard errors are:
rescale(Pc = Pc, std.err = se.Pc, method = "duotrio")

## Can also do
rescale(Pd = c(.6,.7), std.err = c(.2, NA))
psyfun(2, method = "triangle")
psyinv(0.8, method = "twoAFC")
psyderiv(2, method = "duotrio")
pc2pd(0.7, 1/2)
pd2pc(0.3, 1/3)

}
\keyword{models}% at least one, from doc/KEYWORDS

