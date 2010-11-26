`duotrio` <-
  function () 
{
  duotrio <- binomial()
  duotrio$link <- "Link for the duo-trio test"
  duotrio$linkinv <- function(eta) {
    tres <- 1 - pnorm(eta/sqrt(2)) - pnorm(eta/sqrt(6)) + 
      2 * pnorm(eta/sqrt(2)) * pnorm(eta/sqrt(6))
    tres[eta < 0] <- 0.5
    tres
  }
  duotrio$mu.eta <- function(eta) -dnorm(eta/sqrt(2))/sqrt(2) - 
    dnorm(eta/sqrt(6))/sqrt(6) + 2 *
      (dnorm(eta/sqrt(2)) * pnorm(eta/sqrt(6))/sqrt(2) +
       pnorm(eta/sqrt(2)) * dnorm(eta/sqrt(6))/sqrt(6))
  duotrio$linkfun <- function(mu) {
    duotriog <- function(d, p) -p + duotrio$linkinv(d)
    tres <- mu
    for(i in 1:length(mu)) {
      if(mu[i] <= 0.5) 
        tres[i] <- 0
      else if(mu[i] > 1-1e-9)
        tres[i] <- Inf
      else 
        tres[i] <- uniroot(duotriog, c(0, 14), p = mu[i])$root
    }
    tres
  }
  duotrio
}

`threeAFC` <-
  function () 
{
  threeAFC <- binomial()
  threeAFC$link <- "Link for the 3-AFC test"
  threeAFC$linkinv <- function(eta) {
    threeAFCg <- function(x, d) dnorm(x - d) * pnorm(x)^2
    tres <- eta
    for (i in 1:length(eta)) {
      if (eta[i] > 0) 
        tres[i] <- integrate(threeAFCg, -Inf, Inf, d = eta[i])$value
      else tres[i] <- 1/3
    }
    tres
  }
  threeAFC$mu.eta <- function(eta) {
    threeAFCgd <- function(x, d) (x - d) * dnorm(x - d) * 
      pnorm(x)^2
    tres <- eta
    for (i in 1:length(eta)) tres[i] <-
      integrate(threeAFCgd, -Inf, Inf, d = eta[i])$value
    tres
  }
  threeAFC$linkfun <- function(mu) {
    threeAFCg2 <- function(d, p) -p + threeAFC$linkinv(d)
    tres <- mu
    for (i in 1:length(mu)) {
      if (mu[i] > 1/3) 
        res <- uniroot(threeAFCg2, c(0, 10), p = mu[i])$root
      if (mu[i] <= 1/3) 
        res <- 0
      tres[i] <- res
    }
    tres
  }
  threeAFC
}

`triangle` <-
  function () 
{
  triangle <- binomial()
    triangle$link <- "Link for the triangle test"
  triangle$linkinv <- function(eta) {
    triangleg <- function(x, d) 2 * dnorm(x) *
      (pnorm(-x * sqrt(3) + d * sqrt(2/3)) +
       pnorm(-x * sqrt(3) - d * sqrt(2/3)))
    tres <- eta
    for (i in 1:length(eta)) {
      if (eta[i] > 0) 
        tres[i] <- integrate(triangleg, 0, Inf, d = eta[i])$value
      else tres[i] <- 1/3
    }
    tres
  }
  triangle$mu.eta <- function(eta) {
    trianglegd <- function(x, d) 2 * sqrt(2/3) * dnorm(x) * 
      (dnorm(-x * sqrt(3) + d * sqrt(2/3)) -
       dnorm(-x * sqrt(3) - d * sqrt(2/3)))
    tres <- eta
    for (i in 1:length(eta)) tres[i] <-
      integrate(trianglegd, 0, Inf, d = eta[i])$value
    tres
  }
  triangle$linkfun <- function(mu) {
    triangleg2 <- function(d, p) -p + triangle$linkinv(d)
    tres <- mu
    for (i in 1:length(mu)) {
      if (mu[i] > 1/3) 
        tres[i] <- uniroot(triangleg2, c(0, 10), p = mu[i])$root
      if (mu[i] <= 1/3) 
        tres[i] <- 0
    }
    tres
  }
  triangle
}

twoAFC <-
  function () 
{
  twoAFC <- binomial()
  twoAFC$link <- "Link for the 2-AFC test"
  twoAFC$linkfun <- function(mu) {
    tres <- mu
    for (i in 1:length(mu)) {
      if (mu[i] > 0.5) 
        tres[i] <- sqrt(2) * qnorm(mu[i])
      if (mu[i] <= 0.5) 
        tres[i] <- 0
    }
    tres
  }
  twoAFC$linkinv <- function(eta) {
    tres <- pnorm(eta/sqrt(2))
    tres[eta < 0] <- 0.5
    tres
  }
  twoAFC$mu.eta <- function(eta) dnorm(eta/sqrt(2))/sqrt(2)
  twoAFC
}
