`findcr` <-
function (sample.size, alpha = .05, p0 = .5) 
{
  ## Find the critical value of a one-tailed binomial test.
  m <- match.call(expand.dots=FALSE)
  m[[1]] <- as.name("list")
  eval.parent(m) # evaluate the *list* of arguments
  ss <- sample.size 
  if(ss != trunc(ss) | ss <= 0)
    stop("'sample.size' has to be a positive integer")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  if(p0 <= 0 | p0 >= 1)
    stop("'p0' has to be between zero and one")
  ## Core function:
  i <- 0
  while (1 - pbinom(i, ss, p0) > alpha) i <- i + 1
  i + 1
}

