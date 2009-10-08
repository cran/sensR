\name{betabin}
\alias{betabin}
\alias{summary.betabin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Beta-Binomial model for over-dispersed binomial data}
\description{
  Fits the beta binomial model and the chance corrected beta binomial
  model to binomial data.
}
\usage{
betabin(data, start = c(.5,.5), method = c("mu-gamma", "alpha-beta"),
        vcov = TRUE, corrected = TRUE, pGuess = 1/2,
        gradTol = 1e-4, ...)

\method{summary}{betabin}(object, alpha=.05, ...)

}
\arguments{
  \item{object}{an object of class "betabin", i.e. the result of
    \code{betabin()}.}
  \item{alpha}{the allowed type I error for confidence intervals}
  \item{data}{matrix or data.frame with two columns; first column
    contains the number of success and the second the total number of
    cases. The number of rows should correspond to the number of
    observations.}
  \item{start}{starting values to be used in the optimization}
  \item{vcov}{logical, should the variance-covariance matrix of the
    parameters be computed?}
  \item{method}{the desired representation. Optimization is performed in
    the "mu-gamma" parameterization and subsequently transformed to the
    "alpha-beta" scale if needed.}
  \item{corrected}{should the chance corrected or the standard beta
    binomial model be estimated?}
  \item{pGuess}{the guessing probability in the test protocol, i.e. the
    point at which the probability scale is truncated in the chance
    corrected model. The argument is further used in the likelihood
    ratio test of association reported by the summary method for both
    models.}
  \item{gradTol}{a warning is issued if max|gradient| < gradTol, where
    'gradient' is the gradient at the values at which the optimizer
    terminates. This is not used as a termination or convergence
    criterion during model fitting.}
  \item{\dots}{\code{betabin}: The only recognized (hidden) argument is
    \code{doFit} (boolean) which by default is \code{TRUE}. When
    \code{FALSE} \code{betabin} returns an environment which facilitates
    examination of the likelihood surface via the (hidden) functions
    \code{getParBB} and \code{setParBB}.
    Not used in \code{summary.betabin}.}
}
\details{
The following additional methods are implemented objects of class
\code{betabin}:
\code{print}, \code{vcov}, \code{logLik} and \code{coef}.

The bounded optimization is performed with the "L-BFGS-B" optimizer in
\code{\link{optim}}.

Observe that standard errors and confidence intervals are often not
appropriate for parameter close to the boundaries of the parameter
space. The likelihood ratio tests reported by the summary method are
more accurate.

Variance-covariance matrix (and standard errors) is based on the inverted
Hessian at the optimum. The Hessian is obtained with the \code{hessian}
function from the numDeriv package.

The summary method provides a likelihood ratio test of over-dispersion
on one degree of freedom and a likelihood ratio test of association
(i.e. of the hypothesis if "no difference") on two degrees of
freedom (chi-square tests). Since the gamma parameter is tested on the
boundary of the
parameter space, the correct degree of freedom for the first test is
probably 1/2 rather than one, or somewhere in between, and the latter
test is probably also on less than two degrees of freedom. Research is
needed to determine the appropriate no. degrees of freedom to use in
each case. The choices used here are believed to be conservative, so the
stated p-values are probably a little too large.

The parameters of the standard (i.e. corrected = FALSE) beta-binomial
model refers to the mean (i.e. probability) and dispersion on the scale
of the observations, i.e. on the scale where we talk of a probability of
a correct answer. The parameters of the chance corrected (i.e. corrected
= TRUE) beta-binomial model refers to the mean and dispersion on the
scale of "proportion of discriminators" in the population of
interest. The mean parameter (mu) is therefore restricted to the
interval from zero to one in both models, but the has different
interpretations.

The log-likelihood of the standard beta-binomial model is
\deqn{\ell(\alpha, \beta; x, n) = - N \log Beta(\alpha, \beta) +
  \sum_{j=1}^N \log Beta(\alpha + x_j, \beta - x_j + n_j)}

and the log-likelihood of the chance corrected beta-binomial model is
\deqn{\ell(\alpha, \beta; x, n) = - N \log Beta(\alpha, \beta) +
  \sum_{j=1}^N \log \left\{ \sum_{i=1}^{x_j} {{x_j} \choose i}
  (1-p_g)^{n_j-x_j+i} p_g^{x_j-i}
  Beta(\alpha + i, n_j - x_j + \beta) \right\}}{%
  \ell(\alpha, \beta; x, n) = - N \log Beta(\alpha, \beta) +
  \sum_{j=1}^N \log ( \sum_{i=1}^x_j {x_j \choose i}
  (1-p_g)^{n_j-x_j+i} p_g^{x_j-i}
  Beta(\alpha + i, n_j - x_j + \beta) )}

where \eqn{\mu = \alpha/(\alpha + \beta)},
\eqn{\gamma = 1/(\alpha + \beta + 1)}, \eqn{Beta} is the Beta
function, cf. \code{\link{beta}},
\eqn{N} is the number of independent binomial observations, i.e. the
number of rows in \code{data}, and \eqn{p_g} is the guessing
probability, \code{pGuess}.

The gradient at the optimum is evaluated with \code{gradient} from the
numDeriv package.

}
\value{
  An object of class \code{betabin} with elements
  \item{coefficients}{named vector of coefficients}
  \item{vcov}{variance-covariance matrix of the parameter estimates if
    \code{vcov = TRUE}}
  \item{data}{the data supplied to the function}
  \item{call}{the matched call}
  \item{logLik}{the value of the log-likelihood at the MLEs}
  \item{method}{the method used for the fit}
  \item{convergence}{0 indicates convergence. For other error messages,
    see \code{?optim}.}
  \item{message}{possible error message - see \code{?optim} for
    details}
  \item{counts}{the number of iterations used in the optimization - see
    \code{?optim} for details}
  \item{corrected}{Is the chance corrected model estimated?}
  \item{logLikNull}{log-likelihood of the binomial model with
    prop = pGuess}
  \item{logLikMu}{log-likelihood of a binomial model with
    prop = sum(x)/sum(n)}
}
\references{Brockhoff, P.B. (2003). The statistical power of
  replications in difference tests.
  Food Quality and Preference, 14, pp. 405--417.
}

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
(bb <- betabin(dat, corrected = FALSE, method = "mu-gamma"))
summary(bb)
vcov(bb)
logLik(bb)
AIC(bb)
coef(bb)

}
\keyword{models}% at least one, from doc/KEYWORDS
