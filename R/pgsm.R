pgsm <-
function(q, weight, rateparam) {
  numcomp <- length(weight)
  mixcomp <- matrix(NA, nrow = length(q), ncol = numcomp)
  for (i in 1:numcomp) mixcomp[,i] <- pgamma(q, shape = i, rate = rateparam)
  cumprob <- mixcomp%*%weight
  return(cumprob)
}

