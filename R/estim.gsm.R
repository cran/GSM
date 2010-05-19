estim.gsm <-
function(y, J, G = 100, M = 600, a, b, alpha) {
  N <- length(y)
  y.grid <- seq(min(y)*.66, max(y)*1.5, length = G)

  lbl <- matrix(NA, nrow = N, ncol = M)
  wgt <- matrix(NA, nrow = M, ncol = J)

  x <- rep(1, N)
  w <- rep(1 / J, J)
  theta <- rep(NA, M)
  logfdens <- matrix(NA, J, G)
  fdens <- matrix(NA, M, G)

  for (m in 1:M) {
    for (nn in 1:N) {
      temp <- log(seq(a + sum(x) - x[nn], a + sum(x) - x[nn] + J - 1))
      pi <- (1:J - 1)*log(y[nn]) - lgamma(1:J) + cumsum(temp) - (1:J)*log(b + sum(y))
      pi <- (w*exp(pi)) / sum(w*exp(pi))
      x[nn] <- sample(1:J, 1, prob = pi)
      lbl[nn,m] <- x[nn]
    }

    x.counts <- table(factor(x, levels = 1:J))
    w <- rdirichlet(1, x.counts + alpha)
    theta[m] <- rgamma(1, sum(x) + a, rate = sum(y) + b)
    for (j in 1:J) {
      logfdens[j,] <- log(w[j]) + (j - 1)*log(y.grid) - theta[m]*y.grid - lgamma(j) + j*log(theta[m])
    }
    fdens[m,] <- apply(exp(logfdens), 2, sum)
    wgt[m,] <- w
    if (m / 100 == round(m / 100)) print(paste("simulation", m, "/", M))
  }
  result <- list(fdens = fdens, theta = theta, weight = wgt, data = y)
  class(result) <- "gsm"
  return(result)
}

