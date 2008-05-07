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
discrimPwr(delta, sample.size, alpha = .05,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{The underlying sensory difference (non-negative)}
  \item{sample.size}{The sample size (a positive integer)}
  \item{alpha}{The type I level of the test (must be between zero and one)}
  \item{method}{The discrimination test protocol. Four allowed values:
    \code{"twoAFC"}, \code{"threeAFC"}, \code{"duotrio"},
    \code{"triangle"}}  
}
\details{
  The function uses one of the dedicated binomial families.
}
\value{
  The power.
}
\references{Brockhoff, P.B. and Christensen, R.H.B (2008). Thurstonian
  models for sensory discrimination tests as generalized linear models. 
  Manuscript for Food Quality and Preference.}
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
}
\keyword{models}% at least one, from doc/KEYWORDS
