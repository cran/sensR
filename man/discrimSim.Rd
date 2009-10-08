\name{discrimSim}
\alias{discrimSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulates replicated difference tests}
\description{
  Simulates the outcome of \code{sample.size} replicated sensory
  difference tests (for any one of four protocols: 2-AFC, 3-AFC, duotrio
  and triangle tests)
  for a given d-prime value and a given overdispersion (default 0).
}
\usage{
discrimSim(sample.size, replicates, delta, sd.indiv = 0,
           method = c("duotrio", "halfprobit", "probit", "triangle",
           "twoAFC", "threeAFC"))
}
\arguments{
  \item{sample.size}{The sample size - number of subjects}
  \item{replicates}{Number of replications per subject}
  \item{delta}{The delta (d-prime) value}
  \item{method}{The discrimination protocol}
  \item{sd.indiv}{the individual variability in d-prime values. A value
    of 0 (default) is default and corresponds to complete independence}
}
\details{
  The function uses the four thurstonian binomial families.
}
\value{
  A vector of length \code{sample.size} with the number of correct answers for
  each subject.
}
\references{Brockhoff, P.B. and Christensen, R.H.B. (2010). Thurstonian
  models for sensory discrimination tests as generalized linear models.
  Food Quality and Preference, 21, pp. 330-338.}
\author{Rune Haubo B Christensen and Per Bruun Brockhoff}

\seealso{ \code{\link{triangle}}, \code{\link{twoAFC}},
  \code{\link{threeAFC}}, \code{\link{duotrio}},
  \code{\link{discrimPwr}}, \code{\link{discrim}},
  \code{\link{AnotA}}, \code{\link{discrimSS}},
  \code{\link{samediff}}, \code{\link{findcr}}}

\examples{

## Running simulations:
discrimSim(sample.size = 10, replicates = 3, delta = 2,
           method = "triangle", sd.indiv = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
