\name{betabin}
\alias{betabin}
\alias{summary.betabin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Beta-Binomial model to overdispersed data}
\description{
  Fits the beta binomial model to data. 
}
\usage{
betabin(data, start = c(.5,.5), method = c("mu-gamma", "alpha-beta"),
        vcov = TRUE, ...)

\method{summary}{betabin}(object, alpha=.05, ...)

}
\arguments{
  \item{object}{an object of class "betabin", ie. the result of
    \code{betabin()}.}
  \item{alpha}{the allowed type I error for confidence intervals}
  \item{data}{matrix or data.frame with two columns; first column
    contains the number of success and the second the total numnber of
    cases. The number of rows should correspond to the number of
    observations.} 
  \item{start}{starting values to be used in the optimization}
  \item{vcov}{logical, should the variance-covariance matrix of the
    parameters be computed?}
  \item{method}{The desired representation. Note that while the
    "mu-gamma" is often the most natural and easiest to interpret, it
    can cause convergence problems when overdispersion is close to
    non-existence.} 
  \item{\dots}{additional arguments passed to \code{optim} in
    \code{betabin}. Not used in \code{summary.betabin}.}
}
\details{
The following additional methods are implemented objects of class
\code{betabin}: 
\code{print}, \code{vcov}, \code{logLik} and \code{coef}. 
}
\value{
  An object of class \code{betabin} with elements
  \item{coef}{named vector of coefficients}
  \item{vcov}{variance-covariance matrix of the parameter estimates}
  \item{data}{a named vector with the data supplied to the function}
  \item{call}{the matched call}
  \item{logLik}{the value of the log-likelihood at the MLEs}
  \item{method}{the method used for the fit}
  \item{convergence}{0 indicates convergence. For other error messages,
    see \code{?optim}.}
  \item{message}{possible error messsage - see \code{?optim} for
    details}
  \item{counts}{the number of iterations used in the optimization - see
    \code{?optim} for details}
}
\references{Brockhoff, P.B. (2003). The statistical power of
  replications in difference tests.}
\author{Rune Haubo B Christensen}

\seealso{ \code{\link{triangle}}, \code{\link{twoAFC}},
  \code{\link{threeAFC}}, \code{\link{duotrio}},
  \code{\link{discrimPwr}}, \code{\link{discrimSim}},
  \code{\link{discrimSS}}, \code{\link{samediff}},
  \code{\link{AnotA}}, \code{\link{findcr}}
}

\examples{
## Create data:
x <- c(3,2,6,8,3,4,6,0,9,9,0,2,1,2,8,9,5,7)
n <- c(10,9,8,9,8,6,9,10,10,10,9,9,10,10,10,10,9,10)
dat <- data.frame(x, n)

(bb <- betabin(dat, method = "mu-gamma"))
summary(bb)
vcov(bb)
logLik(bb)
AIC(bb)
coef(bb)

}
\keyword{models}% at least one, from doc/KEYWORDS
