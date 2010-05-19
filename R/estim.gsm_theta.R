estim.gsm_theta <-
function(y, J, G = 100, M = 600, a, b, alpha) {
  N <- length(y)
  y.grid <- seq(min(y)*.66, max(y)*1.5, length = G)

  lbl <- matrix(NA, nrow = N, ncol = M)
  wgt <- matrix(NA, nrow = M, ncol = J)

  x <- rep(1, N)
  w <- rep(1 / J, J)
  theta <- rep(J / max(y), M + 1)
  logfdens <- matrix(NA, J, G)
  fdens <- matrix(NA, M, G)

  for (m in 1:M) {
    for (nn in 1:N) {
      pi <- (1:J - 1)*log(y[nn]) - theta[m]*y[nn] - lgamma(1:J) + (1:J)*log(theta[m])
      pi <- (w*exp(pi)) / sum(w*exp(pi))
      x[nn] <- sample(1:J, 1, prob = pi)
      lbl[nn,m] <- x[nn]
    }

    x.counts <- table(factor(x, levels=1:J))
    w <- rdirichlet(1, x.counts + alpha)
    theta[m+1] <- rgamma(1, sum(x) + a, rate = sum(y) + b)
    for (j in 1:J) {
      logfdens[j,] <- log(w[j]) + (j - 1)*log(y.grid) - theta[m + 1]*y.grid - lgamma(j) + j*log(theta[m + 1])
    }
    fdens[m,] <- apply(exp(logfdens), 2, sum)
    wgt[m,] <- w
    if (m / 100 == round(m / 100)) print(paste("simulation", m, "/", M))
  }
  result <- list(fdens = fdens, theta = theta, weight = wgt, data = y)
  class(result) <- "gsm"
  return(result)
}

