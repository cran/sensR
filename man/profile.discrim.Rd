\name{profile.discrim}
\alias{profile.discrim}
\alias{plot.profile.discrim}
\alias{confint.discrim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Profile likelihood methods for discrim objects.}
\description{
  Computes the (normalized or relative) profile likelihood for the
  parameters of a discrimination test, plots the normalized profile
  likelihood and computes profile likelihood confidence intervals.
}
\usage{
\method{profile}{discrim}(fitted, min = 0, max = 3, numpts = 50, ...)

\method{plot}{profile.discrim}(x, level = c(0.99, 0.95), fig = TRUE,
            method = "natural", n = 500, ...)

\method{confint}{discrim}(object, parm, level = 0.95, ...) 
}
\arguments{
  \item{fitted}{a \code{discrim} object.}
  \item{x}{a \code{profile.discrim} object.}
  \item{object}{a \code{discrim} object.}
  \item{parm}{currently not used.}
  \item{min}{the minimum delta for which to do the profiling. By default
    set to 0, which for numerical stability is change internally to 1e-4.}
  \item{max}{the maximum delta beyond the MLE for which to do the profiling.}
  \item{numpts}{control parameter: At how many points should the profile
    likelihood be evaluated?}
  \item{method}{the type of spline to be used in approximating the
    profile likelhood curve (trace)---se \code{\link{spline}} for
    details.}
  \item{n}{the number of spline interpolations to use in plotting the
    profile likelihood curve (trace).}
  \item{level}{for \code{plot}: At which levels to include horizontal lines to indicate
    confidence levels in plots of the normalized profile
    likelihoods. For \code{confint}: at which level to compute the
    confidence interval.}
  \item{fig}{logical: Should the normalized profile likelihoods be plotted?}
  \item{\dots}{For \code{plot}: additional arguments to \code{plot}. For
    \code{confint}: additional arguments to \code{confint.glm} in
    package MASS. For \code{profile}: additional arguments to
    \code{\link{glm}}.}
}
%% \details{
%%   Note that the variance-covariance matrix does not always exist in
%%   contrast to the profile likelihood. \code{profile = FALSE} may
%%   therefore cause confidence intervals etc. to be \code{NA}.
%% }
\value{
  For \code{profile}:
  An object of class \code{"profile.discrim", "data.frame"}---a
  \code{data.frame} with two columns giving
  the value of the parameter and the corresponding value of the profile
  likelihood.

  For \code{plot}:
  An object of class \code{"nProfile.discrim", "data.frame"}---the
  \code{data.frame} from the \code{profile}-object with an extra
  columns containing the
  normalized profile likelihood.

  For \code{confint}:
  
  A 2x2 matrix with columns named \code{"lower", "upper"} giving the
  lower and upper (1 - \code{alpha})\% confidence interval for the
  parameters named in the rows.
}
\references{Brockhoff, P.B. and Christensen R.H.B.(2008). Thurstonian
  models for sensory 
  discrimination tests as generalized linear models.
  Manuscript for Food Quality and Preference.}
\author{Rune Haubo B Christensen and Per Bruun Brockhoff}

%% \seealso{
%%   \code{\link{summary.discrim}}
%% }

\examples{
## 7 success out of 10 samples in a duo-trio experiment:
dd <- discrim(7, 10, "duotrio")
plot(profile(dd))
confint(dd)
points(confint(dd), rep(.1465, 2), pch = 3)

}
\keyword{models}% at least one, from doc/KEYWORDS
