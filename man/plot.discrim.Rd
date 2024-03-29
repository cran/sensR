\name{plot.discrim}
\alias{plot.discrim}
\title{ Plot function for discrim objects}
\description{
  This function plots the latent distributions of sensory intensity
  corresponding to the items or products tested in the discrimination
  test.   
}
\usage{
\method{plot}{discrim}(x, main = TRUE, length = 1000, ...)
}
\arguments{
  \item{x}{The \code{discrim} object whose latent distributions are to
    be plotted}
  \item{main}{include an automatically generated title on the plot?
    Default is \code{TRUE}}
  \item{length}{the length of the vectors to be plotted. Longer vectors
    gives more smooth curves.}
  \item{\dots}{additional arguments to \code{plot} and \code{lines}}
}
\value{
   The function produces a plot and does not return any value.
}
\author{Rune Haubo B Christensen}
\examples{
## Generate discrim objects to be plotted:
fm1 <- discrim(10, 15, method = "threeAFC")
fm2 <- discrim(10, 15, method = "triangle")
old <- par(mfrow=c(2,1)) ## Split plotting window in two
## Plot the distributions of sensory intensity for the two objects
## and increase the line width
plot(fm1, lwd=2) 
plot(fm2, lwd=2)
par(old)
}
\keyword{hplot}
