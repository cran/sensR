`findcr` <-
    function (sample.size, alpha = .05, p0 = .5, pd0 = 0,
              type = c("difference", "similarity"))
{
  ## Find the critical value of a one-tailed binomial test.
    type <- match.arg(type)
    ss <- sample.size
    if(ss != trunc(ss) | ss <= 0)
        stop("'sample.size' has to be a positive integer")
    if(alpha <= 0 | alpha >= 1)
        stop("'alpha' has to be between zero and one")
    if(p0 <= 0 | p0 >= 1)
        stop("'p0' has to be between zero and one")
    if(pd0 < 0 | pd0 > 1)
        stop("'pd0' has to be between zero and one")
    ## Core function:
    i <- 0
    if(type == "difference") {
        while (1 - pbinom(i, ss, pd0 + p0*(1-pd0)) > alpha) i <- i + 1
        i + 1
    }
    else {
        while(pbinom(i, ss, pd0 + p0*(1-pd0)) < alpha) i <- i + 1
        i - 1
    }
}

