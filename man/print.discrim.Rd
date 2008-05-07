\name{print.discrim}
\alias{print.discrim}
\title{Print method for objects of class "discrim"}
\description{
  The function will print results from a discrimination test stored in
  an object of class \code{"discrim"} made with one of the functions
  \code{\link{AnotA}}, \code{\link{discrim}}.
}
\usage{
\method{print}{discrim}(x, digits = getOption("digits"), alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"discrim"}}
  \item{digits}{number of digits in resulting table of results}
  \item{alpha}{the allowed type I level of the Wald based confidence
    interval}
  \item{...}{currently not used}
}
\value{
  A table of results of the descrimination test is printed including
  estimates, standard errors, lower confidence limit, upper confidence
  limit and a p-value for a test of the null-hypothesis of the parameter
  being zero. 
}
\references{Brockhoff, P.B. and Christensen, R.H.B.(2008). Thurstonian
  models for sensory discrimination tests as generalized linear models. 
  Manuscript for Food Quality and Preference.}
\author{Rune Haubo B Christensen}
\seealso{\code{\link{AnotA}}, \code{\link{discrim}},
  \code{\link{discrimPwr}}, \code{\link{discrimSim}},
  \code{\link{discrimSS}}, \code{\link{findcr}}} 

\examples{
## Print results of an A-Not A test:
AnotA(10, 20, 3, 20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
\keyword{models}% __ONLY ONE__ keyword per line
