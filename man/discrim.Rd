\name{discrim}
\alias{discrim}
% \alias{discrimr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensory discrimination analysis}
\description{
  Computes the d-prime, its standard error and P-value for the
  hypothesis of no sensory 
  difference for any  of four methods: 2-AFC, 3-AFC, duotrio and
  triangle tests. 
}
\usage{
discrim(success, total,
        method = c("duotrio", "threeAFC", "twoAFC", "triangle"), ...)  

%% discrimr(formula, data, weights, start, subset, na.action,
%%          contrasts = NULL, method = c("duotrio", "probit", "threeAFC",
%%            "triangle", "twoAFC", "logit"), Hess = TRUE, ...)
}
\arguments{
  \item{success}{the number of correct answers}
  \item{total}{the total number of answers (sample size)}
%%%   \item{formula}{an object of class "formula" (or one that can be coerced
%%%     to that class): a symbolic description of the model to be fitted.}
%%%   \item{data}{an optional data frame, list or environment (or object coercible
%%%     by as.data.frame to a data frame) containing the variables in the
%%%     model. If not found in data, the variables are taken from
%%%     environment(formula)}
%%%   \item{start}{starting values to be used in the optimization}
%%%   \item{subset}{an optional vector specifying a subset of observations to be
%%%     used in the fitting process.} 
%%%   \item{weights}{an optional vector of weights to be used in the fitting
%%%     process. Should be NULL or a numeric vector.}
%%%   \item{na.action}{a function which indicates what should happen when the
%%%     data contain NAs. The default is set by the na.action setting of
%%%     options, and is na.fail if that is unset. The "factory-fresh" default
%%%     is na.omit. Another possible value is NULL, no action. Value na.exclude
%%%     can be useful.} 
%%%   \item{contrasts}{an optional list of specifying the contrasts to use.}
%%%   \item{hess}{logical, should the hessian be computed? (requires for
%%%     standard errors etc. to be computed by summary.discrimr)}
  \item{method}{the discrimination protocol. Four allowed values:
    "twoAFC", "threeAFC", "duotrio", "triangle"} 
  \item{\dots}{additional arguments passed to \code{\link{glm}}}
}
\details{
  The function uses the one of the dedicated binomial families and the
  base function \code{\link{glm}} %% (discrim) or \code{optim} (discrimr). 
}
\value{
  An object of class \code{discrim} with elements
  \item{coef}{ named vector of coefficients (d-prime and tau)}
  \item{res.glm}{the glm-object from the fitting process}
  \item{vcov}{variance-covariance matrix of the coefficients}
  \item{se}{named vector with standard error of the coefficients
    (standard error of d-prime}
  \item{data}{a named vector with the data supplied to the function}
  \item{p.value}{one-sided p-value from Fisher's exact test
    (\code{fisher.test})} 
  \item{test}{a string with the name of the test (\code{A-Not A}) for
    the print method}
  \item{call}{the matched call}
}
\references{Brockhoff, P.B. and Christensen, R.H.B.(2008). Thurstonian
  models for sensory discrimination tests as generalized linear models.
  Manuscript for Food Quality and Preference.}
\author{Rune Haubo B Christensen and Per Bruun Brockhoff}

\seealso{ \code{\link{triangle}}, \code{\link{twoAFC}},
  \code{\link{threeAFC}}, \code{\link{duotrio}},
  \code{\link{discrimPwr}}, \code{\link{discrimSim}},
  \code{\link{discrimSS}}, \code{\link{samediff}},
  \code{\link{AnotA}}, \code{\link{findcr}}
}

\examples{
## Running the simple discrimination test analyses:
discrim(10, 15, "twoAFC")
discrim(10, 15, "threeAFC")
discrim(10, 15, "duotrio")
discrim(10, 15, "triangle")
}
\keyword{models}% at least one, from doc/KEYWORDS
