\name{clm2twoAC}
\alias{clm2twoAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract 2-AC coefficient table from a cumulative link model}
\description{
  The Thurstonian model for the 2-AC protocol can be formulated as a
  cumulative link model (see the references). This function extracts the
  2-AC model parameter 
  estimates, standard errors, z-value and p-values from a cumulative
  link (mixed) model fitted with \code{\link[ordinal]{clm}} or
  \code{\link[ordinal]{clmm}} from package \code{ordinal}.
}
\usage{

clm2twoAC(object, ...)

}
\arguments{
  \item{object}{ a \code{\link[ordinal]{clm}} or
    \code{\link[ordinal]{clmm}} object
    }
  \item{\dots}{not currently used.}
}
\seealso{
  \code{\link{twoAC}}, \code{\link{twoACpwr}}
}
\value{
  A \code{data.frame} with the coefficient table. The two first rows
  contain the estimates of \code{tau} and \code{d.prime} while the
  remaining rows contain optional regression variables for
  \code{d.prime}.
}
\references{
  Christensen R.H.B., Lee H-S and Brockhoff P.B. (2011). Estimation of
  the Thurstonian model for the 2-AC protocol. Submitted to Food
  Quality and Preference. 
}

\author{Rune Haubo B Christensen}

\examples{

## Example of a simple 2-AC model. First the conventional way:
twoAC(c(2, 2, 6))
## The using a cumulative link model (clm from package ordinal): 
library(ordinal)
response <- gl(3,1)
fit.clm <- clm(response ~ 1, weights = c(2, 2, 6), link = "probit")
clm2twoAC(fit.clm)
## Alternatively we could get estimates and standard errors "by hand":
tab <- coef(summary(fit.clm))
theta <- tab[,1]
(tau <- (theta[2] - theta[1])/sqrt(2))
(d.prime <- (-theta[2] - theta[1])/sqrt(2))
VCOV <- vcov(fit.clm)
(se.tau <- sqrt((VCOV[1,1] + VCOV[2,2] - 2*VCOV[2,1])/2))
(se.d.prime <- sqrt((VCOV[1,1] + VCOV[2,2] + 2*VCOV[2,1])/2))

## Extended example with a regression model for d.prime
## (see the referenced paper for details):
n.women <- c(2, 2, 6)*10
n.men <- c(1, 2, 7)*10
wt <- c(n.women, n.men)
response <- gl(3,1, length = 6)
gender <- gl(2, 3, labels = c("women", "men"))
fm2 <- clm(response ~ gender, weights = wt, link = "probit")
clm2twoAC(fm2)

}
\keyword{models}% at least one, from doc/KEYWORDS
