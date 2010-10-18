`AnotA` <-
function (x1, n1, x2, n2, ...)
{
  m <- match.call(expand.dots = FALSE)
  m[[1]] <- as.name("c")
  m <- eval.parent(m) # evaluate the *list* of arguments
  data <- m
  for(i in data) {
    if(i != abs(trunc(i)) | i==0)
      stop("Data have to be positive integers")
  }
  if(x1 >= n1)
    stop("'x1' has to be smaller than 'n1'")
  if(x2 >= n2)
    stop("'x2' has to be smaller than 'n2'")
  call <- match.call()
  test <- "A-Not A"
  ## Arrange data:
  xt <-  cbind(c(x1, x2), c(n1 - x1, n2 - x2))
  ## Fit GLM:
  res <- glm(xt ~ gl(2,1), family = binomial(link = probit), ...)
  ## Prepare output:
  b <- coef(summary(res))
  coef <- -b[2,1] # d-prime
  se <- b[2,2]
  vcov <- as.matrix(se^2) # variance-covariance
  p.value <- fisher.test(xt, alternativ="greater")$p.value
  ## Naming:
  names(vcov) <- names(se) <- names(coef) <- "d-prime"
  fit <- list(coefficients = coef, res.glm = res, vcov = vcov, se = se,
              data = data, p.value = p.value, test = test,
              call = call)
  class(fit) <- c("discrim")
  fit
}

## `discrimPwr` <-
##   function (delta, sample.size, alpha = 0.05,
##             method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
##             pd0 = 0, type = c("difference", "similarity"))
## {
## ### m <- match.call(expand.dots=FALSE)
##     method <- match.arg(method)
##     type <- match.arg(type)
## ### m[[1]] <- as.name("list")
## ### m$method <- NULL
## ### eval.parent(m) # evaluate the *list* of arguments
##     ss <- sample.size
##     ## Control arguments:
##     if(ss != trunc(ss) | ss <= 0)
##         stop("'sample.size' has to be a positive integer")
##     if(delta < 0) stop("'delta' has to be non-negative")
##     if(alpha <= 0 | alpha >= 1)
##         stop("'alpha' has to be between zero and one")
##     if(pd0 < 0 | pd0 > 1)
##         stop("'pd0' has to be between zero and one")
##     ## Find the prob corresponding to delta:
##     prob <- switch(method,
##                    duotrio = duotrio(),
##                    triangle = triangle(),
##                    twoAFC = twoAFC(),
##                    threeAFC = threeAFC() )
##     ## prob under the null hypothesis:
##     Pguess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
##     ## critical value in a one-tailed binomial test:
##     xcr <- findcr(ss, alpha, Pguess, type = type, pd0)
##     ## Compute power of the test:
##     if(type == "difference")
##         power <- 1 - pbinom(xcr - 1, ss, prob$linkinv(delta))
##     else
##         power <- pbinom(xcr, ss, prob$linkinv(delta))
##     power
## }
## 

discrimPwr <-
  function(delta, sample.size, alpha = 0.05,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
           pd0 = 0, type = c("difference", "similarity"))
{
  ## match and test arguments:
  method <- match.arg(method)
  type <- match.arg(type)
  ss <- sample.size
  if(ss != trunc(ss) | ss <= 0)
    stop("'sample.size' has to be a positive integer")
  if(delta < 0) stop("'delta' has to be non-negative")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  if(pd0 < 0 | pd0 > 1)
    stop("'pd0' has to be between zero and one")
  ## get Pc from delta:
  Pc <- psyfun(delta, method = method)
  Pguess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  ## critical value in one-tailed binomial test:
  xcr <- findcr(ss, alpha, Pguess, type = type, pd0)
  ## compute power of the test from critical value:
  if(type == "difference") {
    xcr <- delimit(xcr, low = 1, up = ss + 1)
    power <- 1 - pbinom(q = xcr - 1, size = ss, prob = Pc)
  }
  else if(type == "similarity") {
    xcr <- delimit(xcr, low = 0, up = ss)
    power <- pbinom(q = xcr, size = ss, prob = Pc)
  }
  else ## should never happen
    stop("'type' not recognized")
  return(power)
}

`discrimSS` <-
  function (delta, power, alpha = 0.05,
            method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
            pd0 = 0, type = c("difference", "similarity"), start = 1)
{
    method <- match.arg(method)
    type <- match.arg(type)
    if(delta < 0) stop("'delta' has to be non-negative")
    if(alpha <= 0 | alpha >= 1)
        stop("'alpha' has to be between zero and one")
    if(power <= 0 | power >= 1)
        stop("'power' has to be between zero and one")
    if(pd0 < 0 | pd0 > 1)
        stop("'pd0' has to be between zero and one")
    i <- start
    while (discrimPwr(delta = delta, sample.size = i,
                      alpha = alpha, method = method,
                      pd0 = pd0, type = type) < power)
        i <- i + 1
    i
}

`discrim` <-
function (success, total,
          method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
          pd0 = 0, type = c("difference", "similarity"), ...)
{
    m <- match.call(expand.dots=FALSE)
    method <- match.arg(method)
    type <- match.arg(type)
    m[[1]] <- as.name("list")
    m$method <- m$type <- NULL
    m <- eval.parent(m) # evaluate the *list* of arguments
    x <- m$success;  n <- m$total
    call <- match.call()
    if(x != trunc(x) | x<=0)
        stop("'success' has to be a positive integer")
    if(n != trunc(n) | n<=0)
        stop("'total' has to be a positive integer")
    if(x >= n)
        stop("'total' has to be larger than 'success'")
    if(pd0 < 0 | pd0 > 1)
        stop("'pd0' has to be between zero and one")
    p <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
    ## Compute p-value:
    p.value <-
        if(type == "difference")
            1 - pbinom(x - 1, n, pd0 + p * (1 - pd0))
        else
            pbinom(x, n, pd0 + p * (1 - pd0))
    if(x/n > p) {
##         stop("'succes'/'total' has to be larger than ",
##              ifelse(p < 1/2, "1/3", "1/2") , " for the ",
##              method, " test")
    ## Create glm-family object
        fam <- switch(method,
                      duotrio = duotrio(),
                      threeAFC = threeAFC(),
                      twoAFC = twoAFC(),
                      triangle = triangle() )
        ## Compute d-prime:
        xt <- matrix(c(x, n - x), ncol = 2)
        etastart <- fam$linkfun(xt[1, 1]/sum(xt))
        res <- glm(xt ~ 1, family = fam, etastart = etastart,
                   control = glm.control(epsilon = 1e-05, maxit = 50), ...)
        ## Prepare output:
        s.res <- summary(res)
        coef <- coef(res)
        vcov <- s.res$cov.unscaled
        se <- s.res$coef[,2]
    }
    else {
        coef <- 0
        vcov <-  se <- NA
        res <- NULL
    }
    data <- c(success = x, total = n)
    ## Naming:
    names(coef) <- names(se) <- names(vcov) <- "d-prime"
    ## Output'ing and class'ing:
    fit <- list(coefficients = coef, res.glm = res, vcov = vcov, se = se,
                data = data, p.value = p.value, test = method,
                call = call)
    class(fit) <- "discrim"
    fit
}

`discrimSim` <-
  function(sample.size, replicates, delta, sd.indiv = 0,
           method = c("duotrio", "halfprobit", "probit", "triangle",
           "twoAFC", "threeAFC"))
{
  method <- match.arg(method)
###   m <- match.call(expand.dots=FALSE)
###   m[[1]] <- as.name("list")
###   m$method <- NULL
###   eval.parent(m) # evaluate the *list* of arguments
  method <- match.arg(method)
  if(sample.size != trunc(sample.size) | sample.size <= 0)
    stop("'sample.size' has to be a positive integer")
  if(replicates != trunc(replicates) | replicates < 0)
    stop("'replicates' has to be a non-negativ integer")
  if(delta < 0) stop("'delta' has to be non-negative")
  if(sd.indiv < 0) stop("'sd.indiv' has to be non-negative")
  D <- rnorm(sample.size, 0, sd.indiv) # individual deviations
  q <- delta + D # individual delta's
  ## Compute binomial probabilities:
  ## Note that q[q < 0] <- 0 is not neccesary; this is handled in the
  ## link functions.
  prob <- switch(method,
                 probit = pnorm(q),
                 halfprobit = .5 + .5 * pnorm(q),
                 duotrio = duotrio()$linkinv(q),
                 triangle = triangle()$linkinv(q),
                 twoAFC = twoAFC()$linkinv(q),
                 treeAFC = threeAFC()$linkinv(q))
  ## Compute no. correct answers for the test:
  n.correct <- rbinom(sample.size, replicates, prob)
  n.correct
}

print.discrim <- function(x, digits = getOption("digits"),
                           alpha=.05, ...) {
  ## Print function for the same-diff method.
  coef <- x$coef; se <- x$se; p.value <- x$p.value
  p <- 1-alpha/2
  lower <- coef - qnorm(p) * se
  lower <- ifelse(lower <= 0, 0, lower)
  upper <- coef + qnorm(p) * se
  mat <- c(coef, se, lower, upper, p.value)
  table <- matrix(mat, nrow = length(coef))
  rownames(table) <- names(coef)
  colnames(table) <- c("Estimate", "Std. Error", "Lower",
                       "Upper", "P-value")
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Results for the", x$test, "test:\n\n")
  print(table, digits = digits)
  invisible(x)
}

plot.discrim <-
  function(x, main = TRUE, length = 1000, ...)
{
  z <- seq(-5, 5, length.out = length)
  y <- dnorm(z)
  y2 <- dnorm(z, mean = x$coef[1])
  main.txt <- ifelse(main,
                     paste("Distribution of sensory intensity for the",
                           x$test, "test"), c("") )
  plot(z, y, type="l", xlab = "Sensory Magnitude",
###       ylab = "Sensory Intensity",
       ylab = "", main = main.txt, las = 1, lty = 2, ...)
  lines(z, y2, col = "red", lty = 1, ...)
  invisible()
}

## summary.discrim <- function() {
##
## }
##
## print.summary.discrim <- function() {
##
## }
## fitted.discrim <- function() {
##
## }

discrimr <- ## Discrim revised
  function(formula, data, weights, start, subset, na.action, contrasts
           = NULL, method = c("duotrio", "probit", "threeAFC",
                     "triangle", "twoAFC", "logit"), Hess = TRUE, ...)
{
  nll <- function(beta, X, y, w) { # negative log-likelihood
    eta <- offset
    eta <- eta + X %*% beta
    p <- link.inv(eta)
    if(all(p > 0 && p < 1))
      -sum(2 * w * (y*log(p/(1 - p)) + log(1-p)))
    else(Inf)
  }
  grd <- function(beta, X, y, w) { # gradient
    eta <- offset
    eta <- eta + X %*% beta
    p <- link.inv(eta)
    if(all(p > 0)) {
      ##      cat(NROW(w), NROW(p), NROW(y), "X = ", dim(X), "\n")
      -2 * t(X) %*% (w * muEta(eta) * (y/p - (1-y)/(1-p)))
      ##  browser()
    }
    else(rep(NA, length(beta)))
  }
  m <- match.call(expand.dots = FALSE)
  m$start <- m$Hess <- m$method <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  if (is.matrix(eval.parent(m$data)))
    m$data <- as.data.frame(data)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  n <- nrow(x)
  cons <- attr(x, "contrasts")
  wt <- model.weights(m)
  if (!length(wt))
    wt <- rep(1, n)
  offset <- model.offset(m)
  if (length(offset) <= 1)
    offset <- rep(0, n)
  y <- model.response(m)
  if (NCOL(y) == 2) {
    n <- y[, 1] + y[, 2]
    y <- ifelse(n == 0, 0, y[, 1]/n)
    wt <- wt * n
  }
  stopifnot(all(wt > 0))
  if(missing(start)) start <- rep(0, ncol(x))
  if(missing(data))
  if(length(start) != ncol(x))
    stop("'start' is not of correct length")
  method <- match.arg(method)
  dn <- dimnames(x)[[2]]
  link.inv <- switch(method,
                     duotrio = duotrio()$linkinv,
                     probit = pnorm,
                     threeAFC = threeAFC()$linkinv,
                     triangle = triangle()$linkinv,
                     twoAFC = twoAFC()$linkinv,
                     logit = plogis
###                     twoAFC = function(eta) {pnorm(eta/sqrt(2))}
                     )
  muEta <- switch(method,
                  duotrio = duotrio()$mu.eta,
                  probit = dnorm,
                  threeAFC = threeAFC()$mu.eta,
                  triangle = triangle()$mu.eta,
                  twoAFC = twoAFC()$mu.eta,
                  logit = dlogis
                  )

  fit <- optim(start, fn=nll, gr=grd, X=x, y=y, w=wt, method="BFGS", #
               hessian = Hess, ...)
  ## Fitted values (probabilities):
  fit$fitted <- p <- link.inv(offset + x %*% fit$par)
  ## Deviance:
  fit$dev <- 2 * wt * (y * log(y/p) + (1 - y) * log((1 - y)/(1 - p)))
  fit$deviance <-  sum(fit$dev)
  fit$resid.dev <- sign(y - p) * sqrt(fit$dev)
  se <- NULL
  if(Hess) {
    fit$vcov <- solve(fit$hessian)
    fit$se <- sqrt(diag(fit$vcov))
    names(fit$se) <- dn
    dimnames(fit$vcov) <- list(dn, dn)
  }
  fit$coef <- fit$par
  fit$data <- data
  fit$test <- method
  fit$call <- match.call()
  names(fit$coef) <- dn
  class(fit) <- "discrimr"
  fit
}

confint.discrim <- function(object, parm, level = 0.95, ...)
{
  MASS:::confint.glm(object = object$res.glm, parm = parm,
                     level = level, ...)
}

profile.discrim <-
  function(fitted, min = 0, max = 3, numpts = 50, ...)
{
  dev <- double(numpts)
  da <- fitted$data
  xt <- cbind(da[1], da[2] - da[1])
  min <- ifelse(min < 1e-4, 1e-4, min)
  delta <- seq(min, fitted$coef + max, length = numpts)
  fam <- family(fitted$res.glm)
  for(i in 1:numpts)
    dev[i] <- glm(xt ~ -1 + offset(delta[i]), family = fam, ...)$deviance
  pll <- -dev/2
  Data <- data.frame(delta = delta, pll = pll)
  class(Data) <- c("profile.discrim", "data.frame")
  Data
}

plot.profile.discrim <-
  function(x, level = c(0.99, 0.95), fig = TRUE, method = "natural",
           n = 500, ...)
{
  max.pll <- max(x$pll)
  x$npl <- exp(x$pll - max.pll)
  lim <- sapply(level, function(x) exp(-qchisq(x, df = 1)/2))
  if (fig == TRUE) {
    npl.spline <- spline(x$delta, x$npl, n = n, method = method)
    plot(npl.spline$x, npl.spline$y, type = "l", las = 1, ylim = c(0, 1),
         xlab = expression(delta), ylab = "Normalized Profile Likelihood",
         main = "", ...)
    abline(h = lim, col = "grey")
  }
  class(x) <- c("nProfile.discrim", "data.frame")
  invisible(x)
}
