predict.gsm <-
function(object, thresh, start = 1, ...) {
  numsim <- dim(object$fdens)[1]
  if (numsim < start)
  	stop(paste("number of draws to use smaller than those available (", numsim, ")", sep = ""))
  thetasim <- object$theta[start:numsim]
  weightsim <- object$weight[start:numsim, ]
  runs <- dim(weightsim)[1]
  prob <- vector(length = runs)
  for (i in 1:runs) prob[i] <- 1 - pgsm(thresh, weightsim[i, ], thetasim[i])
  prob
}

