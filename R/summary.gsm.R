summary.gsm <-
function(object, plot = FALSE, start = 1, ...) {
  numsim <- dim(object$fdens)[1]
  if (numsim < start)
  	stop(paste("number of draws to use smaller than those available (", numsim, ")", sep = ""))
  y <- object$data
  thetasim <- object$theta[start:numsim]
  weightsim <- object$weight[start:numsim, ]
  if (numsim == start) {
    out <- list(theta = summary(thetasim), `weights posterior means` = weightsim)
  }
  else {
    out <- list(theta = summary(thetasim), `weights posterior means` = colMeans(weightsim))
  }
  if (plot) {
    J <- dim(object$weight)[2]
    if (numsim == start) {
      barplot(weightsim, names.arg = 1:J, xlab = "Mixture component", ylab = "Estimated posterior mean", main = "Mixture weights posterior means")
    }
    else {
      barplot(colMeans(weightsim), names.arg = 1:J, xlab = "Mixture component", ylab = "Estimated posterior mean", main = "Mixture weights posterior means")
    }
  }
  out
}

