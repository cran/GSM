estim.gsm <- function(y, J, G = 100, M = 600, a, b, alpha) {
  rdirich <- function (n, alpha) {
      l <- length(alpha)
      x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
      sm <- x %*% rep(1, l)
      x/as.vector(sm)
  }

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
    w <- rdirich(1, x.counts + alpha)
    theta[m] <- rgamma(1, sum(x) + a, rate = sum(y) + b)
    for (j in 1:J) {
      logfdens[j,] <- log(w[j]) + (j - 1)*log(y.grid) - theta[m]*y.grid - lgamma(j) + j*log(theta[m])
    }
    fdens[m,] <- apply(exp(logfdens), 2, sum)
    wgt[m,] <- w
    if (m / 100 == round(m / 100)) print(paste("simulation", m, "/", M))
  }
  out <- new("gsm", fdens = fdens, theta = theta, weight = wgt, data = y)
  return(out)
}

estim.gsm_theta <- function(y, J, G = 100, M = 600, a, b, alpha) {
  rdirich <- function (n, alpha) {
      l <- length(alpha)
      x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
      sm <- x %*% rep(1, l)
      x/as.vector(sm)
  }

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
    w <- rdirich(1, x.counts + alpha)
    theta[m+1] <- rgamma(1, sum(x) + a, rate = sum(y) + b)
    for (j in 1:J) {
      logfdens[j,] <- log(w[j]) + (j - 1)*log(y.grid) - theta[m + 1]*y.grid - lgamma(j) + j*log(theta[m + 1])
    }
    fdens[m,] <- apply(exp(logfdens), 2, sum)
    wgt[m,] <- w
    if (m / 100 == round(m / 100)) print(paste("simulation", m, "/", M))
  }
  out <- new("gsm", fdens = fdens, theta = theta, weight = wgt, data = y)
  return(out)
}

allcurves.q <- function(post, perc) {
  n <- dim(post)[2]
  temp <- rep(NA,n)
  for (i in 1:n) temp[i] <- quantile(post[ , i] ,perc)
  return(temp)
}

dgsm <- function(x, weight, rateparam) {
  numcomp <- length(weight)
  mixcomp <- matrix(NA, nrow = length(x), ncol = numcomp)
  for (i in 1:numcomp) mixcomp[ , i] <- dgamma(x, shape = i, rate = rateparam)
  dens <- mixcomp%*%weight
  return(dens)
}

pgsm <- function(q, weight, rateparam) {
  numcomp <- length(weight)
  mixcomp <- matrix(NA, nrow = length(q), ncol = numcomp)
  for (i in 1:numcomp) mixcomp[,i] <- pgamma(q, shape = i, rate = rateparam)
  cumprob <- mixcomp%*%weight
  return(cumprob)
}

rgsm <- function(n, weight, rateparam) {
  J <- length(weight)
  rmixg <- vector(length = n)
  tmp.rmixg <- matrix(NA, nrow = J, ncol = n)
  tmp.lbl <- rmultinom(n, 1, weight)
  for (i in 1:J) {
  	tmp.rmixg[i,] <- rgamma(n, shape = i, rate = rateparam)
  	rmixg <- rmixg + tmp.lbl[i,]*tmp.rmixg[i,]
  }
  return(rmixg)
}
