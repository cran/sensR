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
discrimSS(delta, power, alpha,
          method = c("duotrio", "threeAFC", "twoAFC", "triangle"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{the underlying sensory difference (larger than zero)}
  \item{power}{the wanted power (between zero and one)}
  \item{alpha}{the type 1 level of the test (between zero and one)}
  \item{method}{the discrimination protocol. Four allowed values:
    "twoAFC", "threeAFC", "duotrio", "triangle"}  
}
\details{
  The function uses one of the dedicated binomial families and the
  \code{discrimPwr} function 
}
\value{
  The sample size
}
\references{Brockhoff, P.B. and Christensen, R.H.B.(2008). Thurstonian
  models for sensory discrimination tests as generalized linear models. 
  Manuscript for Food Quality and Preference.} 
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

}
\keyword{models}% at least one, from doc/KEYWORDS
