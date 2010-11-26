\name{SDT}
\alias{SDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Signal Detection Theory Computation of d-prime}
\description{
  The function computes d-prime for any 2 x J table where J >= 2 for the
  "yes--no" or "A-Not A" experiment using the Signal Detection Theory
  (SDT) algorithm to compute J-1 d-prime's. The algorithm is also called
  the "empirical probit transform". The function also provides the "logit"
  counterpart.
}
\usage{
SDT(tab, method = c("probit", "logit"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A 2 x J table with true class relation in rows (only two
    true classes) and the J-class response in columns}
  \item{method}{should the empirical probit or logit transform be computed?}
}
%% \details{
%%   ~~ If necessary, more details than the description above ~~
%% }
\value{
  A (J-1) x 3 matrix. The first two columns contains the z-transform of
  the Hit rate and the False Alarm rate respectively---ready to plot
  along with the empirical ROC curve. The third column contains the
  estimated d-primes.
}
\references{MacMillan , A. N. and Creelman, C. D (2005) Detection Theory
A User's Guide. Lawrence Elbaum Associates, Inc. 2nd edition.}
\author{Rune Haubo B Christensen}
%% \note{ ~~further notes~~ 
%% 
%%  ~Make other sections like Warning with \section{Warning }{....} ~
%% }
%% \seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
### Design table:
## 8  "yes"-responses to no-samples
## 1  "yes"-responses to yes-samples
## 17 "no"-response to no-samples
## 24 "no"-responses to yes-samples
## Note that response-class is columnwise and true-class is rowwise. 
(mat <- matrix(c(8, 17, 1, 24), 2, byrow = TRUE))
SDT(mat, "logit")
SDT(mat, "probit")
## compare to AnotA():
AnotA(8, 25, 1, 25)
## Multi-response-class example (odor example from MacMillan and
## Creelman, 2005)
(odor <- matrix(c(112, 112, 72, 53, 22, 4, 7, 38, 50, 117, 101, 62), 2,
               byrow = TRUE))
obj <- SDT(odor)
ROC(obj[3,3])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
