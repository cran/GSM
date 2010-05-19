allcurves.q <-
function(post, perc) {
  n <- dim(post)[2]
  temp <- rep(NA,n)
  for (i in 1:n) temp[i] <- quantile(post[ , i] ,perc)
  return(temp)
}

