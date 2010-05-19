dgsm <-
function(x, weight, rateparam) {
  numcomp <- length(weight)
  mixcomp <- matrix(NA, nrow = length(x), ncol = numcomp)
  for (i in 1:numcomp) mixcomp[ , i] <- dgamma(x, shape = i, rate = rateparam)
  dens <- mixcomp%*%weight
  return(dens)
}

