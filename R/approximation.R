
# calculates the asymptotic p-value using the method of Robbins and Pitman (1949)
.pAsymptotic <- function(x, lams, eps, acc) {
# x: quantile
# lams: vector of eigenvalues
# eps: accuracy
# acc: accuracy for removing small eigenvalues

  lams <- .weed(lams = lams, accuracy = acc)
  lams <- sort(lams, decreasing=TRUE)
  m <- length(lams)
  if (m == 0)
    p <- NA
  else {
    bet <- min(lams)

    # get an upper bound to the number of iterations needed
    Q2 <- qnorm(eps)^2
    maxiter <- trunc(0.5 * (x/bet + Q2 + sqrt(2*Q2*x/bet + Q2*Q2) - m))

    # starting values
    d <- numeric(maxiter)
    c <- numeric(maxiter + 1)
    c[1] <- prod(sqrt(bet / lams))
    restc <- 1 - c[1]
    chi <- pchisq(x / bet, df = m, lower.tail = FALSE)
    partialsum <- c[1] * chi
    dbase <- (1 - bet / lams)
    ready <- FALSE
    ix <- 1

    # iterate
    while (!ready) {
      d[ix] <- 0.5 * sum(dbase^ix)
      c[ix+1] <- mean(c[1:ix] * d[ix:1])
      if (restc > 100 * .Machine$double.neg.eps) {
        restc <- restc - c[ix+1]
      } else {
        restc <- c[ix+1] * lams[1] / lams[m]
      }
      chi <- pchisq(x / bet, df = m + 2 * ix, lower.tail = FALSE)
      partialsum <- partialsum + c[ix+1] * chi
      error <- restc * pchisq(x / bet, df = m + 2 * ix + 2, lower.tail = TRUE)
      ready <- (error < eps) || (ix == maxiter) || (error / partialsum < 10^-4)
      ix <- ix + 1
    }
    p <- partialsum + restc
    if (p < eps) { p <- 0 }
  }
  p
}

################################################################################

# Removes extremely small eigenvalues
.weed <- function(lams, accuracy) {
# lams: vector of eigenvalues
  if (missing(accuracy)) {
    thresh <- 1/50
  } else {
    thresh <- 1/accuracy
  }
  lams <- -sort(-lams)
  m <- length(lams)
  if (lams[1] == 0) {
    lams <- numeric(0)
    m <- 0
  } else {
    while (lams[m] / lams[1] < thresh) {
      q <- m-1
      r <- m-2
      lams[q] <- lams[q] + lams[m]
      while ((r > 0) && (lams[r] < lams[q])) {
        lams[r:q] <- mean(lams[r:q])
        r <- r - 1
      }
      m <- q
    }
    lams <- lams[1:m]
  }
  lams
}

