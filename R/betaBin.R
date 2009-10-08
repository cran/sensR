betabin <-
  function(data, start = c(.5,.5), method = c("mu-gamma", "alpha-beta"),
           vcov = TRUE, corrected = TRUE, pGuess = 1/2,
           gradTol = 1e-4, ...)
{
    m <- match.call(expand.dots = FALSE)
    call <- match.call()
    m$method <- NULL
    m[[1]] <- as.name("list")
    m <- eval.parent(m)
    doFit <- TRUE # A little trick:
    if("doFit" %in% names(m$...)) doFit <- (m$...)$doFit
    if(is.data.frame(m$data)) m$data <- as.matrix(m$data)
    if(!is.matrix(m$data))
        stop("'data' is not a matrix or data.frame")
    if(NCOL(m$data) != 2 || NROW(m$data) < 3)
        stop("'data' should have 2 columns and > 3 rows")
    if(pGuess <= 0 || pGuess >= 1)
        stop("pGuess has to be in the open interval (0, 1)")
    method <- match.arg(method)
    if(any(start < 1e-3) || any(start > 1- 1e-3))
        stop("start has to be in the open interval (0, 1)")
    name <-
        if(method == "mu-gamma") c("mu", "gamma")
        else c("alpha", "beta")
    bbRho <- bbEnvir(parent.frame(), X=as.matrix(m$data),
                     corrected = corrected, pGuess = pGuess,
                     start = start)
    if(!doFit) return(bbRho)
    fit <- optim(getParBB(bbRho), fn = function(par) setParBB(bbRho, par),
                 method = "L-BFGS-B", hessian = FALSE,
                 lower = bbRho$lbounds, upper = bbRho$ubounds,
                 control = list(parscale = c(.01, .01)))
    if(all(fit$par > bbRho$lbounds) && all(fit$par < bbRho$ubounds)) {
        grad <- grad(function(par) setParBB(bbRho, par), fit$par)
        if(max(abs(grad)) > gradTol)
            warning(sprintf("Optimizer terminated with max|gradient|: %e",
                            max(abs(grad))), call. = FALSE) }
    else
        warning("Parameters at boundary occured", call. = FALSE)
    coef <- fit$par
    if(method == "alpha-beta") {
        coef[1] <- fit$par[1] * (1 - fit$par[2]) / fit$par[2]
        coef[2] <- (1 - fit$par[2]) * (1 - fit$par[1]) / fit$par[2]
    }
    names(coef) <- name
    res <- list(coefficients = coef,
                convergence = fit$convergence,
                message = fit$message, counts = fit$counts,
                call = match.call(), data = m$data,
                method = method, corrected = corrected)
    res$logLik <- -fit$value + bbRho$Factor
    res$logLikNull <-
        sum(dbinom(bbRho$x, bbRho$n, prob = pGuess, log = TRUE))
    res$logLikMu <-
        with(bbRho, sum(dbinom(x, n, prob = sum(x)/sum(n),
                               log = TRUE)))
    if(vcov) {
        res$vcov <- matrix(NA, 2, 2)
        names(res$vcov) <- list(name, name)
        if(all(fit$par > bbRho$lbounds) && all(fit$par < bbRho$ubounds)) {
            if(method == "alpha-beta") bbRho$method <- "alpha-beta"
            res$vcov <-
                solve(hessian(function(par) setParBB(bbRho, par), coef))
        }
    }
    class(res) <- c("betabin")
    res
}

bbEnvir <- function(parent, X, corrected, pGuess, start)
{
    rho <- new.env(parent = parent)
    rho$corrected <- corrected
    rho$method <- "mu-gamma"
    rho$X <- X
    rho$n <- X[,2]
    rho$x <- X[,1]
    rho$N <- nrow(X)
    rho$start <- rho$par <- start
    rho$lbounds <- c(1e-6, 1e-6)
    if(corrected) rho$lbounds <- c(1e-6, 1e-3)
    rho$ubounds <- c(1 - 1e-6, 1 - 1e-6)
##     rho$const <- choose(x[1], 0:x[1]) * (1 - pGuess)^(0:x[1]) *
##         pGuess^(-(0:x[1]))
    rho$nllAux <- function(x, a, b)
        log(sum(choose(x[1], 0:x[1]) * (1 - pGuess)^(0:x[1]) *
                pGuess^(-(0:x[1])) *
                beta(a + 0:x[1], b + x[2] - x[1])))
    ## The term choose(x[1], 0:x[1]) does not have to be computed
    ## every time.
    rho$Factor <- sum(lchoose(rho$n, rho$x))
    if(corrected)
        rho$Factor <- rho$Factor +
            sum((rho$n-rho$x) * log(1-pGuess)) + sum(rho$x * log(pGuess))
    rho
}

getParBB <- function(rho) rho$par
setParBB <- function(rho, par)
{
    if(!missing(par))
        rho$par <- par
    with(rho, {
        if(method == "mu-gamma") {
            a <- par[1] * (1 - par[2]) / par[2]
            b <- (1 - par[2]) * (1 - par[1]) / par[2] }
        else {
            a <- par[1]
            b <- par[2] }
        if(corrected)
            N * lbeta(a, b) - sum(apply(X, 1, nllAux, a=a, b=b))
        else
            N * lbeta(a, b) - sum(lbeta(a + x, b - x + n))
    })
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
    if(is.null(object$vcov))
        stop("summary needs vcov in object")
    object$se <- sqrt(diag(object$vcov))
    names(object$se) <- names(object$coef)
    p <- 1-alpha/2
    object$lower <- object$coef - qnorm(p) * object$se
    object$upper <- object$coef + qnorm(p) * object$se
    table <- cbind(object$coef, object$se, object$lower, object$upper)
    rownames(table) <- names(object$coef)
    colnames(table) <- c("Estimate", "Std. Error", "Lower", "Upper")
    object$LR.OD <- 2 * with(object, logLik - logLikMu)
    object$p.value.OD <-
        pchisq(object$LR.OD, df = 1, lower.tail = FALSE)
    object$LR.null <- 2 * with(object, logLik - logLikNull)
    object$p.value.null <-
        pchisq(object$LR.null, df = 2, lower.tail = FALSE)
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
  cat("\nlog-likelihood: ", x$logLik, "\n")
  cat("LR-test of over-dispersion, G^2:", x$LR.OD, "df:", 1,
      "p-value:", x$p.value.OD, "\n")
  cat("LR-test of association, G^2:", x$LR.null, "df:", 2,
      "p-value:", x$p.value.null, "\n")
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

