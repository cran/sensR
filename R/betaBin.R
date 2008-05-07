nll.betabin.ab <-  function(par, x, n) {
  ## Negative log-likelihood for the beta-binomial model in the
  ## alpha-beta parameterization 
  a <- par[1]; b <- par[2]
  tmp <- lbeta(a,b)
  - sum(lbeta(a + x, b - x + n) - tmp)
}

nll.betabin.mg <- function(par, x, n) {
  ## Negative log-likelihood for the beta-binomial model in the
  ## mu-gamma parameterization
  if(any(par <= 0) || any(par >= 1)) return(1e8)
  mu <- par[1]; g <- par[2]
  a <- (mu/g) * (1 - g)
  b <- (1/g) * (1 - g) * (1 - mu)
  tmp <- lbeta(a, b)
  ## Negative log-likelihood of the Beta-binomial distribution
  nll <- - sum(lbeta(a + x, b - x + n) - tmp)
  ##  cat(nll, tmp, a, b, mu, g, "\n")
  nll
}

nll.ccbetabin.ab <- function(par, x, n, C) {
  if(any(par <= 0) || any(par >= 1)) return(Inf)
  term1 <- function(i, n, x, a, b, C)
    (1 - C)^(n - x + i) * C^(x - i) * beta(a + i, n - x + b)
  a <- par[1]; b <- par[2]
  N <- length(n)
  Term1 <- sapply(1:N, term1, n, x, a, b, C)
  nll <- log(sum(Term1)) - lbeta(a, b)
  cat(nll, a, b, "\n")
  -nll
}


### betabin <- function(x, ...) UseMethod("betabin")
### betabin.default <- function(x, ...) betabin.formula(x, ...)
### betabin.matrix <- function(x, ...) betabin.data.frame(x, ...)
###

nll.bin <- function(par, x, n)
  -sum(x * log(par) + (n - x) * log(1 - par))
       

betabin <-
  function(data, start = c(.5,.5), method = c("mu-gamma", "alpha-beta"),
           vcov = TRUE, ...)
{
### FIXME:
### Formula as a possible first argument as well as numeric vector 
### Corrected betabinomial as an option (separate function)    
  m <- match.call(expand.dots = FALSE)
  call <- match.call()
  m$method <- NULL
  m[[1]] <- as.name("list")
  m <- eval.parent(m)
  method <- match.arg(method)
  x <- m$data[,1]
  n <- m$data[,2]
  ## Optimize log-likelihood:
  if(method == "mu-gamma") {
    fit <- optim(start, fn = nll.betabin.mg, x = x, n = n,
                 method = "L-BFGS-B",
                 lower=c(1e-6, 1e-6), upper=c(1 - 1e-6, 1 - 1e-6),
                 hessian = vcov,
                 control = list(parscale = c(1, .01)))
    name <- c("mu", "gamma")
  }
  else {
    fit <- optim(start, fn = nll.betabin.ab, x = x, n = n,
                 method = "BFGS", hessian = vcov, ...)
    name <- c("alpha", "beta")
  }
  ## Extract Output:
  coef <- fit$par
  names(coef) <- name
  if(vcov) {
    vcov <- solve(fit$hessian)
    dimnames(vcov) <- list(name, name)
  }
  else
    vcov <- NULL
  res <- list(coef = coef, vcov = vcov, convergence = fit$convergence,
              logLik = -fit$value, message = fit$message, counts =
              fit$counts, call = call, data = data, method = method)
  class(res) <- c("betabin")
  res
}

print.betabin <-
  function(x, digits = max(3, getOption("digits") - 3), ... )
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)

}

summary.betabin <-
  function(object, alpha=.05, ...)
{
  object$se <- sqrt(diag(object$vcov))
  names(object$se) <- names(object$coef)
  p <- 1-alpha/2
  object$lower <- object$coef - qnorm(p) * object$se
  object$upper <- object$coef + qnorm(p) * object$se
  table <- cbind(object$coef, object$se, object$lower, object$upper)
  rownames(table) <- names(object$coef)
  colnames(table) <- c("Estimate", "Std. Error", "Lower", "Upper")
  if(object$method == "mu-gamma") 
    object$LR.gamma <- 2 *
      (object$logLik + nll.bin(object$coef[1],
                               object$data[,1], object$data[,2]))
  object$p.value.gamma <- pchisq(object$LR.gamma, df = 1, lower = FALSE)
  object$table <- table
  class(object) <- "summary.betabin"
  object
}

print.summary.betabin <-
  function(x, digits = getOption("digits"), alpha=.05, ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  print(x$table)
  conv <- ifelse(x$convergence, "No", "Yes")
  cat("\nlog-likelihood: ", x$logLik, " Fit converged: ", conv, "\n")
  cat("LR-test of gamma:", x$LR.gamma, "p-value:",
      x$p.value.gamma, "\n")
  invisible(x)
}

vcov.betabin <- function(object, ...) {
  object$vcov
}

logLik.betabin <- function(object, ...) {
  val <- object$logLik
  names(val) <- NULL
  attr(val, "nobs") <- sum(object$data[,2])
  attr(val, "df") <- 2
  class(val) <- "logLik"
  val
}

coef.betabin <- function(object, ...) {
  object$coef
}

## Potential functions:
## profile.betabin <- function() {
## 
## }
## 
## plot.profile.betabin <- function() {
## 
## }
## 
## confint.betabin <- function() {
## 
## }
  
