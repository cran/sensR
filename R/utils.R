rescale <-
  function(Pc, Pd, d.prime, std.err, 
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"))
{
  m <- match.call(expand.dots = FALSE)
  m[[1]] <- as.name("list")
  m <- eval.parent(m) # evaluate the *list* of arguments
  arg <- c("Pc", "Pd", "d.prime")
  isPresent <- sapply(arg, function(arg) !is.null(m[[arg]]))
  if(sum(isPresent) != 1)
    stop("One and only one of Pc, Pd and d.prime should be given")
  method <- match.arg(method)
  Pguess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  par <- arg[isPresent]
  if(!is.null(se <- m$std.err)) {
    stopifnot(is.numeric(se) && length(se) == length(m[[par]]))
    stopifnot(all(se[!is.na(se)] > 0))
  }
  if(par == "Pc") {
    Pc <- m[[par]]
    stopifnot(is.numeric(Pc) && all(Pc >= 0) && all(Pc <= 1))
    tooSmall <- Pc < Pguess
    Pc[tooSmall] <- Pguess
    Pd <- pc2pd(Pc, Pguess)
    d.prime <- psyinv(Pc, method = method)
    if(!is.null(se)) {
      se.Pc <- se
      se.Pc[tooSmall] <- NA
      se.Pd <- se.Pc / (1 - Pguess)
      se.d.prime <- se.Pc / psyderiv(d.prime, method = method)
    }
  }
  if(par == "Pd") {
    Pd <- m[[par]]
    stopifnot(is.numeric(Pd) && all(Pd >= 0) && all(Pd <= 1))
    Pc <- pd2pc(Pd, Pguess)
    d.prime <- psyinv(Pc, method = method)
    if(!is.null(se)) {
      se.Pd <- se
      se.Pc <- se.Pd * (1 - Pguess)
      se.d.prime <- se.Pc / psyderiv(d.prime, method = method)
    }
  }
  if(par == "d.prime") {
    stopifnot(is.numeric(d.prime) && all(d.prime >= 0))
    d.prime <- m[[par]]
    Pc <- psyfun(d.prime, method = method)
    Pd <- pc2pd(Pc, Pguess)
    if(!is.null(se)) {
      se.d.prime <- se
      se.Pc <- se * psyderiv(d.prime, method = method)
      se.Pd <- se.Pc / (1 - Pguess)
    } 
  }
  coef <- data.frame(Pc = Pc, Pd = Pd, d.prime = d.prime)
  res <- list(coefficients = coef)
  if(!is.null(se))
    res$std.err <- data.frame(Pc = se.Pc, Pd = se.Pd,
                              d.prime = se.d.prime)
  res$method <- method
  class(res) <- "rescale"
  return(res)
}

print.rescale <- function(x, digits = getOption("digits"), ...)
{
  cat(paste("\nEstimates for the", x$method, "protocol:\n", sep = " "))
  print(coef(x))
  if(!is.null(x$std.err)) {
    cat("\nStandard errors:\n")
    print(x$std.err)
  }
  return(invisible(x))
}

pc2pd <- function(Pc, Pguess)
### Maps Pc to Pd

### arg: Pc: numeric vector; 0 <= Pc <= 1
###      Pguess: the guessing probability; numeric scalar,
###              0 <= Pc <= 1
### res: Pd: numeric vector; 0 <= Pc <= 1
{
  stopifnot(is.numeric(Pguess) && length(Pguess) == 1 &&
            Pguess >= 0 && Pguess <= 1)
  stopifnot(is.numeric(Pc) && all(Pc >= 0) && all(Pc <= 1))
  Pd <- (Pc - Pguess) / (1 - Pguess)
  Pd[Pc <= Pguess] <- 0
  names(Pd) <- names(Pc)
  return(Pd)
}

pd2pc <- function(Pd, Pguess) {
### Maps Pd to Pc
  
### arg: Pd: numeric vector; 0 <= Pc <= 1
###      Pguess: the guessing probability; numeric scalar,
###              0 <= Pc <= 1
### res: Pc: numeric vector; 0 <= Pc <= 1
  stopifnot(is.numeric(Pguess) && length(Pguess) == 1 &&
            Pguess >= 0 && Pguess <= 1)
  stopifnot(is.numeric(Pd) && all(Pd >= 0) && all(Pd <= 1))
  Pc <- Pguess + Pd * (1 - Pguess)
  names(Pc) <- names(Pd)
  return(Pc)
}

psyfun <-
  function(d.prime,
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"))
### Maps d.prime to Pc for sensory discrimination protocols
  
### arg: d.prime: non-negative numeric vector
### res: Pc: numeric vector
{
  method <- match.arg(method)
  stopifnot(all(is.numeric(d.prime)) && all(d.prime >= 0))
  psyFun <- switch(method,
                   duotrio = duotrio()$linkinv,
                   triangle = triangle()$linkinv,
                   twoAFC = twoAFC()$linkinv,
                   threeAFC = threeAFC()$linkinv)
  Pc <- numeric(length(d.prime))
### Extreme cases are not handled well in the links, so we need: 
  OK <- d.prime < Inf
  if(sum(OK) > 0)
    Pc[OK] <- psyFun(d.prime[OK])
  Pc[!OK] <- 1
  names(Pc) <- names(d.prime)
  return(Pc)
}

psyinv <- function(Pc, 
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"))
### Maps Pc to d.prime for sensory discrimination protocols

### arg: Pc: numeric vector; 0 <= Pc <= 1
### res: d.prime: numeric vector
{
  method <- match.arg(method)
  stopifnot(all(is.numeric(Pc)) && all(Pc >= 0) && all(Pc <= 1))
  psyInv <- switch(method,
                   duotrio = duotrio()$linkfun,
                   triangle = triangle()$linkfun,
                   twoAFC = twoAFC()$linkfun,
                   threeAFC = threeAFC()$linkfun)
  d.prime <- numeric(length(Pc))
### Extreme cases are not handled well in the links, so we need: 
  OK <- Pc < 1
  if(sum(OK) > 0)
    d.prime[OK] <- psyInv(Pc[OK])
  d.prime[!OK] <- Inf
  names(d.prime) <- names(Pc)
  return(d.prime)
}

psyderiv <-
  function(d.prime, 
           method = c("duotrio", "threeAFC", "twoAFC", "triangle"))
### Computes the derivative of the psychometric functions at some
### d.prime for sensory discrimination protocols.
  
### arg: d.prime: non-negative numeric vector
### res: Pc: numeric vector
{
  method <- match.arg(method)
  stopifnot(all(is.numeric(d.prime)) && all(d.prime >= 0))
  psyDeriv <- switch(method,
                     duotrio = duotrio()$mu.eta,
                     triangle = triangle()$mu.eta,
                     twoAFC = twoAFC()$mu.eta,
                     threeAFC = threeAFC()$mu.eta)
  Deriv <- numeric(length(d.prime))
### Extreme cases are not handled well in the links, so we need: 
  OK <- d.prime > 0 && d.prime < Inf
  if(sum(OK) > 0)
    Deriv[OK] <- psyDeriv(d.prime[OK])
  Deriv[d.prime == 0] <- NA
  Deriv[d.prime == Inf] <- 0
  names(Deriv) <- names(d.prime)
  return(Deriv)
}

## findcr <-
##   function (sample.size, alpha = .05, p0 = .5, pd0 = 0,
##             type = c("difference", "similarity"))
## {
##   ## Find the critical value of a one-tailed binomial test.
##   type <- match.arg(type)
##   ss <- sample.size
##   if(ss != trunc(ss) | ss <= 0)
##     stop("'sample.size' has to be a positive integer")
##   if(alpha <= 0 | alpha >= 1)
##     stop("'alpha' has to be between zero and one")
##   if(p0 <= 0 | p0 >= 1)
##     stop("'p0' has to be between zero and one")
##   if(pd0 < 0 | pd0 > 1)
##     stop("'pd0' has to be between zero and one")
##   ## Core function:
##   i <- 0
##   if(type == "difference") {
##     while (1 - pbinom(i, ss, pd0 + p0*(1-pd0)) > alpha) i <- i + 1
##     i + 1
##   }
##   else {
##     while(pbinom(i, ss, pd0 + p0*(1-pd0)) < alpha) i <- i + 1
##     i - 1
##   }
## }

test.crit <- function(xcr, sample.size, p.correct = 0.5, alpha = 0.05, type)
### Is xcr the critical value of a one-tailed binomial test?
### Result: boolean

### OBS: there is deliberately no requirement that xcr should be
### positive or less than sample.size.
{  
  if(type %in% c("difference", "greater")) ## alternative is "greater"
    ((1 - pbinom(q = xcr - 1, size = sample.size, prob = p.correct) <= alpha) &&
     (1 - pbinom(q = xcr - 2, size = sample.size, prob = p.correct) > alpha))
  else if(type %in% c("similarity", "less")) ## alternative is "less"
    ((pbinom(q = xcr, size = sample.size, prob = p.correct) <= alpha) &&
     (pbinom(q = xcr + 1, size = sample.size, prob = p.correct) > alpha))
  else
    stop("unknown type argument")
}

findcr <-
  function(sample.size, alpha = 0.05, p0 = 0.5, pd0 = 0,
           type = c("difference", "similarity"))
### Find the critical value of a one-tailed binomial
### test. "difference" means a "greater" alternative hypothesis and
### "similarity" means a "less" alternative hypothesis.

### FIXME: What should this function do/return if the critical value
### is larger than the sample size? Or when it is negative as can
### happen with similarity? Examples:
### (xcr <- findcr(sample.size = 1, p0 = psyfun(1, "twoAFC"))) ## 2
### (xcr <- findcr(sample.size = 1, type = "similarity")) ## -1
### This means that there is no number large/small enough for this
### sample size that could give a significant p-value. Maybe this
### should just be a deliberate feature.
{
  ## match and test arguments:
  type <- match.arg(type)
  ss <- sample.size
### FIXME: Does this test work as intented?
  if(ss != trunc(ss) | ss <= 0)
    stop("'sample.size' has to be a positive integer")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  if(p0 <= 0 | p0 >= 1)
    stop("'p0' has to be between zero and one")
  if(pd0 < 0 | pd0 > 1)
    stop("'pd0' has to be between zero and one")
  ## core function:
  Pc <- pd2pc(pd0, p0)
  if(type == "difference") {
    crdiff <-  function(cr)
      1 - pbinom(q = cr - 1, size = ss, prob = Pc) - alpha
    interval <- c(0, ss + 2) ## deliberately outside allowed range
  }
  else if(type == "similarity") {
    crdiff <- function(cr)
      pbinom(q = cr + 1, size = ss, prob = Pc) - alpha
    interval <- c(-2, ss) ## deliberately outside allowed range
  }
  else ## should never occur
    stop("'type' not recognized") 
  xcr <- round(uniroot(crdiff, interval = interval)$root)
  ## is xcr the critical value?:
  is.crit <- test.crit(xcr = xcr, sample.size = ss, p.correct = Pc,
                       alpha = alpha, type = type)
  if(is.crit) return(xcr)
  ## if uniroot fails, then do a simple search around the vicinity of
  ## the result from uniroot:
  max.iter <- 20 ## avoid infinite loop
  xcr <- delimit(xcr - 10, low = -1)
  i <- 0
  if(type == "difference") {
    while(1 - pbinom(q = xcr + i, size = ss, prob = Pc) > alpha) {
      if(i > max.iter || xcr + i > ss) break 
      i <- i + 1
    }
    xcr <- xcr + i + 1
  }
  if(type == "similarity") {
    while(pbinom(q = xcr + i, size = ss, prob = Pc) < alpha) {
      if(i > max.iter || xcr + i > ss) break
      i <- i + 1
    }
    xcr <- xcr + i - 1
  }
  ## is xcr now the critical value?:
  is.crit <- test.crit(xcr = xcr, sample.size = ss, p.correct = Pc,
                       alpha = alpha, type = type)
  if(is.crit) return(xcr)
  else stop("Failed to find critical value")
}

delimit <- function(x, lower, upper, set.na = FALSE)
### Sets the values of x < lower to lower and values of x > upper to
### upper. If set.na is TRUE values are set to NA. If both lower and
### upper are supplied, the (lower < upper) has to be TRUE.
{
  m <- match.call()
  m[[1]] <- as.name("list")
  m <- eval.parent(m)
  if(!is.null(m$lower) && !is.null(m$upper))
    stopifnot(m$lower < m$upper)
  if(!is.null(m$lower))
    x[x < m$lower] <- ifelse(set.na, NA, m$lower)
  if(!is.null(m$upper))
    x[x > m$upper] <- ifelse(set.na, NA, m$upper)
  return(x)
}
