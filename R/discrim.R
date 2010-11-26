AnotA <-
  function (x1, n1, x2, n2, ...)
### Revise for better computational methods (don't use glm unless it
### is needed to get the standard errors?) and reconsider the output.
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
### FIXME: Do not use Fishers exact test here!
  p.value <- fisher.test(xt, alternativ="greater")$p.value
  ## Naming:
  names(vcov) <- names(se) <- names(coef) <- "d-prime"
  fit <- list(coefficients = coef, res.glm = res, vcov = vcov, se = se,
              data = data, p.value = p.value, test = test,
              call = call)
  class(fit) <- c("anota")
  fit
}

normalPwr <-
  function(pdA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"))
{
  test <- match.arg(test)
  stopifnot(pdA >= 0 && pdA <= 1)
  stopifnot(pd0 >= 0 && pdA <= 1)
  stopifnot(sample.size > 0)
  stopifnot(alpha > 0 && alpha < 1)
  stopifnot(pGuess >= 0 && pGuess < 1)
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  pc0 <- pd2pc(pd0, pGuess)
  pcA <- pd2pc(pdA, pGuess)
  sigma0 <- sqrt(pc0*(1 - pc0)/sample.size)
  sigmaA <- sqrt(pcA*(1 - pcA)/sample.size)
  if(test == "difference") {
    lambda <- (qnorm(1 - alpha)*sigma0 + pc0 - pcA) / sigmaA
    pwr <- pnorm(lambda, lower = FALSE)
  }
  else if(test == "similarity") {
    lambda <- (qnorm(alpha)*sigma0 + pc0 - pcA) / sigmaA
    pwr <- pnorm(lambda, lower = TRUE)
  }
  else
    stop("'test' not recognized")
  return(pwr)
}

discrimPwr <-
  function(pdA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"),
           statistic = c("exact", "normal"))
{
  ## match and test arguments:
  test <- match.arg(test)
  stat <- match.arg(statistic)
  ss <- sample.size
  stopifnot(length(pdA) == 1 && length(pd0) == 1 &&
            length(sample.size) == 1 && length(alpha) == 1 &&
            length(pGuess) == 1)
  if(ss != trunc(ss) || ss <= 0)
    stop("'sample.size' has to be a positive integer")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  stopifnot(length(pd0) == 1 && pd0 >= 0 && pd0 <= 1)
  ## if(pd0 < 0 | pd0 > 1)
  ##   stop("'pd0' has to be between zero and one")
  if(pdA < 0 | pdA > 1)
    stop("'pdA' has to be between zero and one")
  ## get pc from pdA and pGuess:
  pc <- pd2pc(pdA, pGuess)
  if(stat == "normal") {
    pwr <- normalPwr(pdA = pdA, pd0 = pd0, sample.size = ss,
                     alpha = alpha, pGuess = pGuess, test = test)
    return(pwr)
  } ## stat == "exact":
  ## critical value in one-tailed binomial test:
  xcr <- findcr(ss, alpha, pGuess, test = test, pd0)
  ## compute power of the test from critical value:
  if(test == "difference") {
    xcr <- delimit(xcr, low = 1, up = ss + 1)
    power <- 1 - pbinom(q = xcr - 1, size = ss, prob = pc)
  }
  else if(test == "similarity") {
    xcr <- delimit(xcr, low = 0, up = ss)
    power <- pbinom(q = xcr, size = ss, prob = pc)
  }
  else ## should never happen
    stop("'test' not recognized")
  return(power)
}

d.primePwr <-
  function(d.primeA, d.prime0 = 0, sample.size, alpha = 0.05,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
           test = c("difference", "similarity"),
           statistic = c("exact", "normal"))
{
  ## Convenience function that simply modifies some arguments and
  ## calls discrimPwr
  newCall <- call <- match.call()
  method <- match.arg(method)
  stopifnot(length(d.primeA) == 1 && is.numeric(d.primeA) &&
            d.primeA >= 0)
  stopifnot(length(d.prime0) == 1 && is.numeric(d.prime0) &&
            d.prime0 >= 0)
  pdA <- coef(rescale(d.prime = d.primeA, method = method))$pd
  pd0 <- coef(rescale(d.prime = d.prime0, method = method))$pd
  newCall$method <- newCall$d.primeA <- newCall$d.prime0 <- NULL
  newCall$pGuess <-
    ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  newCall$pdA <- pdA
  newCall$pd0 <- pd0
  newCall[[1]] <- as.name("discrimPwr")
  return(eval.parent(newCall))
}

normalSS <-
  function(pdA, pd0 = 0, target.power = 0.9, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"))
{
  test <- match.arg(test)
  stopifnot(pdA >= 0 && pdA <= 1)
  stopifnot(pd0 >= 0 && pdA <= 1)
  stopifnot(target.power > 0 && target.power < 1)
  stopifnot(alpha > 0 && alpha < 1)
  stopifnot(pGuess >= 0 && pGuess < 1)
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  pc0 <- pd2pc(pd0, pGuess)
  pcA <- pd2pc(pdA, pGuess)
  s0 <- sqrt(pc0*(1 - pc0))
  sA <- sqrt(pcA*(1 - pcA))
  ## approximate sample size:
  if(test == "difference") {
    stopifnot(pdA > pd0)
    n <- ((sA * qnorm(1 - target.power) - s0 * qnorm(1 - alpha)) /
          (pc0 - pcA))^2
  }
  else if(test == "similarity") {
    stopifnot(pdA < pd0)
    n <- ((sA * qnorm(target.power) - s0 * qnorm(alpha)) /
          (pc0 - pcA))^2
  }
  return(ceiling(n))
}

discrimSS <-
  function(pdA, pd0 = 0, target.power = 0.90, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"),
           statistic = c("exact", "normal")) 
{
  test <- match.arg(test)
  call <- match.call()
  stat <- match.arg(statistic)
  stopifnot(length(pdA) == 1 && length(pd0) == 1 &&
            length(target.power) == 1 && length(alpha) == 1 && 
            length(pGuess) == 1)
  if(pdA < 0 | pdA > 1)
    stop("'pdA' has to be between zero and one")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  if(target.power <= 0 | target.power >= 1)
    stop("'target.power' has to be between zero and one")
  if(pd0 < 0 | pd0 > 1)
    stop("'pd0' has to be between zero and one")
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  ssN <- normalSS(pdA = pdA, pd0 = pd0, target.power = target.power,
                  alpha = alpha, pGuess = pGuess, test = test)
  ssN <- ssN - 1
  if(stat == "normal")
    return(ssN + 1)
  ## if(stat == "exact"):
  if(ssN <= 50)
    ssTry <- 1
  if(ssN > 1e4) {
    warning(paste("sample size probably > 1e4 and 'exact' option is",
                  "not available\n",
                  "using normal approximation instead"))
    return(ssN)
  }
  if(ssN > 50 && ssN <= 1e4) {
    for(ssTry in floor(ssN * (9:1 / 10))) {
      tryPwr <-
        discrimPwr(pdA = pdA, pd0 = pd0, sample.size = ssTry,
                   alpha = alpha, pGuess = pGuess, test = test)
      if(tryPwr < target.power)
        break
    }
  }
  while (discrimPwr(pdA = pdA, pd0 = pd0, sample.size = ssTry, 
                    alpha = alpha, pGuess = pGuess, 
                    test = test) < target.power)
    ssTry <- ssTry + 1
  return(ssTry)
}

d.primeSS <- 
  function(d.primeA, d.prime0 = 0, target.power = 0.90, alpha = 0.05,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
           test = c("difference", "similarity"),
           statistic = c("exact", "normal")) 
{
  ## Convenience function that simply modifies some arguments and
  ## calls discrimSS
  newCall <- call <- match.call()
  method <- match.arg(method)
  stopifnot(length(d.primeA) == 1 && is.numeric(d.primeA) &&
            d.primeA >= 0)
  stopifnot(length(d.prime0) == 1 && is.numeric(d.prime0) &&
            d.prime0 >= 0)
  pdA <- coef(rescale(d.prime = d.primeA, method = method))$pd
  pd0 <- coef(rescale(d.prime = d.prime0, method = method))$pd
  newCall$method <- newCall$d.primeA <- newCall$d.prime0 <- NULL
  newCall$pGuess <-
    ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  newCall$pdA <- pdA
  newCall$pd0 <- pd0
  newCall[[1]] <- as.name("discrimSS")
  return(eval.parent(newCall))
}


discrim <-
  function(correct, total, pd0 = 0, conf.level = 0.95,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
           statistic = c("exact", "likelihood", "score", "Wald"),
           test = c("difference", "similarity"), ...)
{
  stopifnot(length(correct) == 1 && length(total) == 1 &&
            length(pd0) == 1 && length(conf.level) == 1)
  m <- match.call(expand.dots=FALSE)
  method <- match.arg(method)
  test <- match.arg(test)
  stat <- match.arg(statistic)
  m[[1]] <- as.name("list")
  m$method <- m$statistic <- m$test <- NULL
  m <- eval.parent(m) # evaluate the *list* of arguments
  x <- m$correct;  n <- m$total
  call <- match.call()
### FIXME: is trunc(x) the right way to do this?
  if(x != trunc(x) || x < 0)
    stop("'correct' has to be a non-negative integer")
  if(n != trunc(n) || n <= 0)
    stop("'total' has to be a positive integer")
  if(x > n)
    stop("'correct' cannot be larger than 'total'")
  if(pd0 < 0 || pd0 > 1)
    stop("'pd0' has to be between zero and one")
  Pguess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  pc0 <- pd2pc(pd0, Pguess)
  if(test == "similarity" && pd0 == 0)
    stop("'pd0' should be positive for a similarity test")
  mu <- x/n
  se.mu <- sqrt(mu*(1 - mu)/n)
  ## Draft coefficient table:
  table <- array(NA, dim = c(3, 4))
  rownames(table) <- c("pc", "pd", "d-prime")
  colnames(table) <- c("Estimate", "Std. Error", "Lower", "Upper")
  ## Fill in estimates:
  obj <- rescale(pc = mu, method = method)
  table[,1] <- unlist(obj$coefficients)
  pc.hat <- table[1,1]
  ## Fill in standard errors:
  if(mu < 1 && mu > Pguess) {
    obj <- rescale(pc = mu, std.err = se.mu, method = method)
    table[,2] <- unlist(obj$std.err)
  }
  ## Get p-value, CI and test statistic:
  if(stat == "exact") {
    p.value <- 
      if(test == "difference")
        1 - pbinom(q = x - 1, size = n, prob = pc0)
      else
        pbinom(q = x, size = n, prob = pc0)
    ci <- c(binom.test(x = x, n = n, alternative = "two.sided",
                       conf.level = conf.level)$conf.int)
  }
  if(stat == "likelihood") {
    prof <- profBinom(x, n, nProf = 100)
    ci <- c(confint(prof, level = conf.level))
    logLikMax <- dbinom(x, n, table[1,1], log = TRUE)
    logLikNull <- dbinom(x, n, pc0, log = TRUE)
    Stat <- sign(table[1,1] - pc0) *
      sqrt(2 * (logLikMax - logLikNull))
    p.value <-
      if(test == "difference") pnorm(Stat, lower = FALSE)
      else pnorm(Stat)
  }
  if(stat == "Wald") {
    Stat <- (pc.hat - pc0) / sqrt(pc.hat*(1 - pc.hat)/n)
    p.value <-
      if(test == "difference") pnorm(Stat, lower = FALSE)
      else pnorm(Stat)
    a <- (1 - conf.level)/2
    ci <- mu + se.mu * qnorm(c(a, 1 - a))
  }
  if(stat == "score") {
    ci <- prop.test(x = x, n = n, alternative = "two.sided",
                    conf.level = conf.level)$conf.int
    if(test == "difference")
      score <- prop.test(x = x, n = n, alternative = "greater",
                         p = pc0, correct = FALSE)
    else
      score <- prop.test(x = x, n = n, alternative = "less",
                         p = pc0, correct = FALSE)
    Stat <- score$statistic
    p.value <- score$p.value
  }
  if(sum(is.na(ci)) == 0) {
    ci <- delimit(x = ci, lower = 0, upper = 1)
    intervals <- rescale(pc = ci, method = method)$coefficients
    table[,3] <- unlist(intervals[1,])
    table[,4] <- unlist(intervals[2,])
  }
  res <- list(coefficients = table, p.value = p.value, call = call,
              test = test, method = method, statistic = stat,
              data = c("correct" = x, "total" = n), pd0 = pd0,
              conf.level = conf.level)
  if(stat != "exact")
    res$stat.value <- Stat
  if(stat == "score")
    res$df <- score$parameter
  if(stat == "likelihood")
    res$profile <- prof
  class(res) <- "discrim"
  return(res)
}
    
discrimOld <-
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
    class(fit) <- "anota"
    fit
}

discrimSim <-
  function(sample.size, replicates, d.prime, sd.indiv = 0,
           method = c("duotrio", "halfprobit", "probit", "triangle",
           "twoAFC", "threeAFC"))
{
  method <- match.arg(method)
  if(sample.size != trunc(sample.size) | sample.size <= 0)
    stop("'sample.size' has to be a positive integer")
  if(replicates != trunc(replicates) | replicates < 0)
    stop("'replicates' has to be a non-negativ integer")
  if(d.prime < 0) stop("'d.prime' has to be non-negative")
  if(sd.indiv < 0) stop("'sd.indiv' has to be non-negative")
  ## Individual deviations in d.prime:
  D <- rnorm(n = sample.size, mean = 0, sd = sd.indiv) 
  q <- delimit(d.prime + D, lower = 0) # individual d.prime's
  ## Compute no. correct answers for the test:
  n.correct <- rbinom(n = sample.size, size = replicates,
                      prob = psyfun(q, method = method))
  n.correct
}

print.anota <-
  function(x, digits = getOption("digits"),
                           alpha=.05, ...) {
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

print.discrim <-
  function(x, digits = getOption("digits"), ...)
{
  text1 <- switch(x$statistic,
                  "exact" = "'exact' binomial test.",
                  "likelihood" = "likelihood root statistic.",
                  "Wald" = "Wald statistic.",
                  "score" = "Pearson and score statistics.")
  cat(paste("\nEstimates for the", x$method,
            "discrimination protocol with", x$data[1],
            "correct\nanswers in", 
            x$data[2], "trials. p-value and",
            round(100 * x$conf.level, 3),  
            "percent confidence intervals\nare based on the", text1,
            "\n\n"))
  print(x$coefficients, digits = digits)
  Pguess <- ifelse(x$method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  d.prime0 <- psyinv(pd2pc(x$pd0, Pguess), method = x$method)
  cat(paste("\nResult of", x$test, "test:\n"))
  if(x$statistic == "Wald")
    cat(paste("Wald statistic =", format(x$stat.value, digits),
              "p-value =", format.pval(x$p.value), "\n"))
  if(x$statistic == "likelihood")
    cat(paste("likelihood root statistic =",
              format(x$stat.value, digits), 
              "p-value =", format.pval(x$p.value), "\n"))
  if(x$statistic == "exact")
    cat(paste("'exact' binomial test:",
              "p-value =", format.pval(x$p.value), "\n"))
  if(x$statistic == "score")
    cat(paste("Pearson X-square statistic = ",
              format(x$stat.value, digits), ", df = ", x$df,
              ", p-value = ", format.pval(x$p.value), "\n", sep = ""))
  cat("Alternative hypothesis: ")
  cat(paste("d-prime is",
            ifelse(x$test == "difference", "greater", "less"),
            "than", format(d.prime0, digits), "\n\n"))
  invisible(x)
}

plot.discrim <-
  function(x, main = TRUE, length = 1000, ...)
{
  z <- seq(-5, 5, length.out = length)
  y <- dnorm(z)
  y2 <- dnorm(z, mean = coefficients(x)[1, 3])
  main.txt <- ifelse(main,
                     paste("Distribution of sensory intensity for the",
                           x$method, "test"), c("") )
  plot(z, y, type="l", xlab = "Sensory Magnitude",
       ylab = "", main = main.txt, las = 1, lty = 2, ...)
  lines(z, y2, col = "red", lty = 1, ...)
  invisible()
}

plot.anota <-
  function(x, main = TRUE, length = 1000, ...)
{
  z <- seq(-5, 5, length.out = length)
  y <- dnorm(z)
  y2 <- dnorm(z, mean = coefficients(x))
  main.txt <- ifelse(main,
                     paste("Distribution of sensory intensity for the",
                           x$test, "test"), c("") )
  plot(z, y, type="l", xlab = "Sensory Magnitude",
       ylab = "", main = main.txt, las = 1, lty = 2, ...)
  lines(z, y2, col = "red", lty = 1, ...)
  invisible()
}

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

confint.anota <-
  function(object, parm, level = 0.95, ...)
### get confint from the discrim object.
{
  ci <- MASS:::confint.glm(object$res, level = level)
  rownames(ci) <- c("threshold", "d.prime")
  ci[2,] <- -ci[2,]
  return(ci)
}

confint.discrim <-
  function(object, parm, level = 0.95, ...)
### get confint from the discrim object.
{
  if(level == object$conf.level) 
    obj <- object$coefficients[, 3:4]
  else {
    call <- update(object, conf.level = level, evaluate = FALSE)
    obj <- eval.parent(call)$coefficients[, 3:4]
  }    
  attr(obj, "method") <- object$method
  attr(obj, "conf.level") <- object$conf.level
  attr(obj, "statistic") <- object$statistic
  return(obj)
}

profile.discrim <-
  function(fitted, ...)
{
  if(fitted$statistic == "likelihood")
    prof <- fitted$profile
  else {
    call <- update(fitted, statistic = "likelihood", evaluate = FALSE)
    fitted <- eval.parent(call)
    prof <- fitted$profile
  }
  pg <- ifelse(fitted$method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  prof <- prof[prof$pSeq >= pg, ]
  prof$d.prime <- psyinv(prof$pSeq, method = fitted$method)
  prof$pSeq <- NULL
  attr(prof, "method") <- fitted$method
  class(prof) <- c("profile.discrim", "data.frame")
  return(prof)
}

plot.profile.discrim <-
  function(x, level = c(0.99, 0.95), fig = TRUE, method = "natural",
           n = 1e3, ...)
{
  lim <- sapply(level, function(x) exp(-qchisq(x, df = 1)/2))
  if (fig == TRUE) {
    npl.spline <- spline(x$d.prime, x$Lroot, n = n, method = method)
    plot(npl.spline$x, exp(-npl.spline$y^2/2), type = "l", las = 1, ylim = c(0, 1),
         xlab = "d-prime", ylab = "Relative Likelihood",
         main = "", ...)
    abline(h = lim)
  }
  ## class(x) <- c("nProfile.discrim", "data.frame")
  invisible(x)
}

profBinom <- function(success, total, nProf = 100, ...)
{
    x <- success; n <- total
    if(x != trunc(x) || x < 0)
        stop("'success' has to be a non-negative integer")
    if(n != trunc(n) || n <= 0)
        stop("'total' has to be a positive integer")
    if(x > n)
        stop("'success' > 'total' not allowed")
    pHat <- x/n
    logLik <- dbinom(x, n, pHat, log=TRUE)
    ## pSeq <- seq(from=1e-4, to=1-1e-4, length=nProf)
    pSeq <- seq(from = 1e-8, to = 1-1e-8, length = nProf)
    ## add 'from' and 'to' to the argument list?
    if(!(pHat %in% pSeq))
        pSeq <- sort(c(pHat, pSeq))
    ll <- dbinom(x, n, pSeq, log=TRUE)
    sign <- 2*(pSeq > pHat) -1
    Lroot <- sign * sqrt(2) * sqrt(-ll + logLik)
    prof <- data.frame(pSeq = pSeq, Lroot=Lroot)
    attr(prof, "logLik") <- logLik
    attr(prof, "pHat") <- pHat ## MLE
    class(prof) <- c("profBinom", "data.frame")
    prof
}

confint.profBinom <- function(object, level=0.95) {
    a <- (1-level)/2
    a <- c(a, 1-a)
    cutoff <- qnorm(a)
    sp <- with(object, spline(plogis(pSeq), Lroot))
    ci <- array(0, dim=c(1,2))
    ci[] <- approx(sp$y, sp$x, xout = cutoff)$y
    ci <- qlogis(ci)
    pHat <- attr(object, "pHat")
    if(pHat == 0) ci[1] <- 0
    if(pHat == 1) ci[2] <- 1
    rownames(ci) <- "p"
    colnames(ci) <- c("lower", "upper")
    attr(ci, "level") <- level
    ci
}

plotProf <-
    function(object, method = c("binom", "duotrio", "threeAFC",
                     "twoAFC", "triangle"), log = FALSE,
             relative = TRUE, ...){
    method <- match.arg(method)
    if(method == "binom") {
        sp <- with(object, spline(pSeq, Lroot))
        xlab <- "p"
    } else {
        fam <- switch(method,
                      duotrio = duotrio(),
                      threeAFC = threeAFC(),
                      twoAFC = twoAFC(),
                      triangle = triangle() )
        dSeq <- fam$linkfun(object$pSeq)
        skip <- seq_len(sum(dSeq == 0) -1) ## skip leading zeros
        sp <- spline(dSeq[-skip], object$Lroot[-skip])
        xlab <- expression(delta)
    }
    if(relative){
        y <- -sp$y^2/2
    } else {
        y <- -sp$y^2/2 + attr(object, "logLik")
    }
    if(!log) y <- exp(y)
    plot(sp$x, y, type ="l", xlab=xlab,
         ylab = "Relative likelihood", ...)
### Make correct ylab for various cases.
    if(relative & !log)
        abline(h = c(0.1465, .0362))
### Include level-argument and compute the levels.
    ## points(ci.pHat, rep(.1465, 2))
### should something be returned? x and y perhaps?
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

## discrimPwr <-
##   function(delta, sample.size, alpha = 0.05,
##            method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
##            pd0 = 0, test = c("difference", "similarity"))
## {
##   ## match and test arguments:
##   method <- match.arg(method)
##   test <- match.arg(test)
##   ss <- sample.size
##   if(ss != trunc(ss) | ss <= 0)
##     stop("'sample.size' has to be a positive integer")
##   if(delta < 0) stop("'delta' has to be non-negative")
##   if(alpha <= 0 | alpha >= 1)
##     stop("'alpha' has to be between zero and one")
##   if(pd0 < 0 | pd0 > 1)
##     stop("'pd0' has to be between zero and one")
##   ## get pc from delta:
##   pc <- psyfun(delta, method = method)
##   Pguess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
##   ## critical value in one-tailed binomial test:
##   xcr <- findcr(ss, alpha, Pguess, test = test, pd0)
##   ## compute power of the test from critical value:
##   if(test == "difference") {
##     xcr <- delimit(xcr, low = 1, up = ss + 1)
##     power <- 1 - pbinom(q = xcr - 1, size = ss, prob = pc)
##   }
##   else if(test == "similarity") {
##     xcr <- delimit(xcr, low = 0, up = ss)
##     power <- pbinom(q = xcr, size = ss, prob = pc)
##   }
##   else ## should never happen
##     stop("'test' not recognized")
##   return(power)
## }


### This is not working due to the discreteness of the binomial
### distribution: 
### testSS <- function(target.power, pdA, pd0, sample.size, alpha = 0.05,
###                    pGuess, test) 
### ### Is sample.size the required sample size for a one-tailed binomial test?
### ### Result: boolean
###   ((discrimPwr2(pdA = pdA, pd0 = pd0, sample.size = sample.size,
###                 alpha = alpha, pGuess = pGuess, test = test)
###     >= target.power) &&
###    (discrimPwr2(pdA = pdA, pd0 = pd0, sample.size = sample.size - 1,
###                 alpha = alpha, pGuess = pGuess, test = test)
###     < target.power))


## discrimSSold <-
##   function (delta, power, alpha = 0.05,
##             method = c("duotrio", "threeAFC", "twoAFC", "triangle"),
##             pd0 = 0, type = c("difference", "similarity"), start = 1)
## {
##     method <- match.arg(method)
##     type <- match.arg(type)
##     if(delta < 0) stop("'delta' has to be non-negative")
##     if(alpha <= 0 | alpha >= 1)
##         stop("'alpha' has to be between zero and one")
##     if(power <= 0 | power >= 1)
##         stop("'power' has to be between zero and one")
##     if(pd0 < 0 | pd0 > 1)
##         stop("'pd0' has to be between zero and one")
##     i <- start
##     while (discrimPwr(delta = delta, sample.size = i,
##                       alpha = alpha, method = method,
##                       pd0 = pd0, type = type) < power)
##         i <- i + 1
##     i
## }

