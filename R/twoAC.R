twoAC <- function(data, ...) {
  nll <- function(par) {
    tau <- par[1]
    d <- par[2]
    p1 <- pnorm(-tau, d, sqrt(2))
    p2 <- pnorm(tau, d, sqrt(2)) - p1
    p3 <- 1 - p1 - p2
    p <- c(p1, p2, p3)
    nll <- -sum(data * log(p))
    nll
  }
  if(!isTRUE(all.equal(round(data), data))) ## use round rather than as.integer
    stop("decimal numbers not allowed in 'data' argument")
  if(length(data) != 3)
    stop("'data' argument should be of length 3")
  coef <- matrix(0, 2, 4, dimnames = list(c("tau", "d.prime"),
                            c("Estimate", "Std. Error", "z value",
                              "Pr(>|z|)")))
  ## if(data[3] <= sum(data[1:2])) {
  ##   d <- coef[1,1] <- 0
  ##   tau <- vcov <- NULL
  ## }
  gamma <- cumsum(data / sum(data))
  z <- qnorm(gamma)[-3]
  z <- z * sqrt(2)
  tau <- (z[2] - z[1]) / 2
  d <- -z[1] - tau
  hess <- hessian(nll, x = c(tau, d))
  vcov <- solve(hess)
  coef[,1] <- c(tau, d)
  coef[,2] <- sqrt(diag(vcov))
  coef[,3] <- coef[,1] / coef[,2]
  coef[,4] <- 2*pnorm(abs(coef[, 3]), lower.tail=FALSE)

  fit <- list(tau = tau, d.prime = d, vcov = vcov, coefficients = coef,
              logLik = -nll(c(tau, d)), data = data, call = match.call())
  class(fit) <- "twoAC"
  fit
}

print.twoAC <- function(x, digits = max(3, .Options$digits - 3),
                        ...) {
  cat("\nCoefficients:\n")
  print(x$coefficients, quote = FALSE, digits = digits, ...)
  cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
  invisible(x)
}

profile.twoAC <-
  function(fitted, alpha = 1e-3, nSteps = 1e2, range, ...)
{
  nll.tau <- function(tau, d, data) {
    p1 <- pnorm(-tau, d, sqrt(2))
    p2 <- pnorm(tau, d, sqrt(2)) - p1
    p3 <- 1 - p1 - p2
    p <- c(p1, p2, p3)
    -sum(data * log(p))
  }
  ## get range from Wald CI or over-rule with range argument:
  if(missing(range))
    ## get range from Wald interval and supplied alpha:
    range <- as.vector(confint(fitted, parm = "d.prime",
                               level = 1 - alpha, type = "Wald"))
  else {
    range <- as.numeric(range)
    ## tjeck that d.prime.hat is in the range interval
    if(fitted$d.prime >= max(range) || fitted$d.prime <= min(range))
      stop("d.prime should be in the interval given by range")
  }
  dseq <- seq(from = min(range), to = max(range),
              length.out = nSteps)
  nll <- ## negative profile log-likelihood
    sapply(dseq, function(dd)
           optimize(nll.tau, c(1e-6, 10), d = dd, data = fitted$data)$objective)
  ## get likelihood root statistic:
  sgn <- 2*(dseq < fitted$d.prime) - 1
  Lroot <- sgn * sqrt(2) * sqrt(nll + fitted$logLik)
  res <- data.frame("Lroot" = c(0, Lroot),
                    "d.prime" = c(fitted$d.prime, dseq))
  res <- res[order(res[,1]),]
  if(!all(diff(res[,2]) < 0))
    warning("likelihood is not monotonically decreasing from maximum,\n",
            "  so profile may be unreliable")
  prof <- vector("list", length = 1)
  names(prof) <- "d.prime"
  prof[[1]] <- res
  val <- structure(prof, original.fit = fitted)
  class(val) <- "profile.twoAC"
  val
}

confint.twoAC <-
  function(object, parm, level = 0.95,
           type = c("profile", "Wald"), ...)
{
  type <- match.arg(type)
  if(type == "Wald") {
    cf <- coef(object)[,1]
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(round(100*a, 1), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    ses <- coef(object)[,2][parm]
    ci[] <- cf[parm] + ses %o% fac
  }  else
  if(type == "profile") {
    object <- profile(object, alpha = (1 - level) / 4, ...) 
    ci <- confint(object, level=level, ...)
  }
  return(ci)
}

confint.profile.twoAC <-
  function(object, parm = "d.prime", level = 0.95, ...)
{
  of <- attr(object, "original.fit")
  if(parm != "d.prime")
    stop("Profile likelihood confidence interval only available for 'd.prime'")
  a <- (1-level)/2
  a <- c(a, 1-a)
  pct <- paste(round(100*a, 1), "%")
  ci <- array(NA, dim = c(1, 2),
              dimnames = list(parm, pct))
  cutoff <- qnorm(a)
  pro <- object[[ "d.prime" ]]
  sp <- spline(x = pro[, 2], y = pro[, 1])
  ci[1, ] <- approx(-sp$y, sp$x, xout = cutoff)$y
  ci
}

plot.profile.twoAC <-
  function(x, level = c(0.95, 0.99), Log = FALSE, relative = TRUE,
           fig = TRUE, n = 1e3, ..., ylim = NULL)
{
  ML <- attr(x, "original.fit")$logLik
  lim <- sapply(level, function(x)
                exp(-qchisq(x, df=1)/2) )
  pro <- x[[ "d.prime" ]]
  sp <- spline(x = pro[, 2], y = pro[, 1], n=n)
  sp$y <- -sp$y^2/2
  if(relative && !Log) {
    sp$y <- exp(sp$y)
    ylab <- "Relative profile likelihood"
    dots <- list(...)
    if(missing(ylim))
      ylim <- c(0, 1)
  }
  if(relative && Log) {
    ylab <- "Relative profile log-likelihood"
    lim <- log(lim)
  }
  if(!relative & Log) {
    sp$y <- sp$y + ML
    ylab <- "profile Log-likelihood"
    lim <- ML + log(lim)
  }
  if(!relative & !Log) {
    stop("Not supported: at least one of 'Log' and 'relative' ",
         "have to be TRUE")
    sp$y <- exp(sp$y + ML)
    ylab <- "Profile likelihood"
    lim <- exp(ML + log(lim))
  }
  x[[ "d.prime" ]] <- sp
  if(fig) {
    plot(sp$x, sp$y, type = "l", ylim = ylim,
         xlab = "d.prime", ylab = ylab, ...)
    abline(h = lim)
  }
  attr(x, "limits") <- lim
  invisible(x)
}

  


