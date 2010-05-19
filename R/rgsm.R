rgsm <-
function(n, weight, rateparam) {
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

