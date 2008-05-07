\name{print.AUC}
\alias{print.AUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print AUC objects}
\description{
  The function prints the estimated AUC of a AUC object as well as
  confidence intervals, they are part of the AUC object.
}
\usage{
\method{print}{AUC}(x, digits = getOption("digits"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The AUC object to be printed}
  \item{digits}{the number of digits desired in printing}
  \item{\dots}{currently not used}
}
%% \details{
%%   ~~ If necessary, more details than the description above ~~
%% }
%% \value{
%%   The function does not return any value
%%   ~Describe the value returned
%%   If it is a LIST, use
%%   \item{comp1 }{Description of 'comp1'}
%%   \item{comp2 }{Description of 'comp2'}
%%   ...
%% }
%% \references{ ~put references to the literature/web site here ~ }
\author{Rune Haubo B Christensen}
%% \note{ ~~further notes~~ 
%% 
%%  ~Make other sections like Warning with \section{Warning }{....} ~
%% }
%% \seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## Obtain d-prime and standard error for an arbitrary test:
fm1 <- AnotA(8, 25, 1, 25)
## Print AUC with and without confidence limits:
AUC(fm1$coef, , fm1$se)
AUC(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
