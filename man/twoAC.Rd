\name{twoAC}
\alias{twoAC}
\alias{print.twoAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{2-AC discrimination protocol}
\description{
  Provides estimates and standard errors of d-prime and tau for the two
  alternative (2-AC) discrimination protocol. The 2-AC protocol is
  equivalent to a 2-AFC protocol with a "no-difference" option. 
}
\usage{

twoAC(data, ...)

}
\arguments{
  \item{data}{
    a numeric vector of length 3 with the number of observations in the
    three response categories. If the third element is larger than the
    first element, then the estimate of d-prime will be positive.
  }
  \item{\dots}{
    not currently used
  }
}
\details{
  \code{\link[=confint.twoAC]{confint}},
  \code{\link[=profile.twoAC]{profile}} and
  \code{print} methods are implemented for \code{twoAC} objects.
}
\value{
  An object of class \code{twoAC} with elements
  \item{d.prime}{
    the estimate of d-prime
  }
  \item{tau}{
    the estimate of tau
  }
  \item{coefficients}{
    named vector of coefficients; d-prime and tau
  }
  \item{vcov}{
    variance-covariance matrix of the parameter estimates
  }
  \item{data}{
    the data supplied to the function
  }
  \item{call}{
    the matched call
  }
  \item{logLik}{
    the value of the log-likelihood at the maximum likelihood estimates
  }
}
\references{
  Christensen, R.H.B. and Brockhoff, P.B. (2010). Estimation of the
  Thurstonian model for the 2-AC protocol.
  Manuscript for Food Quality and Preference.
}

\author{Rune Haubo B Christensen}

\examples{

fit <- twoAC(c(2,2,6))
fit

}
\keyword{models}% at least one, from doc/KEYWORDS
