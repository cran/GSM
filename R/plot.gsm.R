plot.gsm <-
function(x, ndens = 5, xlim = c(min(y), max(y)), ylim = c(0, max(x$fdens)), xlab = "x", ylab = "density", nbin = 10, histogram = FALSE, bands = FALSE, confid = .95, start = 1, ...) {
  numsim <- dim(x$fdens)[1]
  if (numsim < start)
  	stop(paste("number of draws to use smaller than those available (", numsim, ")", sep = ""))
  rpl <- ifelse(ndens > (numsim - start + 1), TRUE, FALSE)
  y <- x$data
  fdens <- x$fdens[start:numsim, ]
  G <- dim(fdens)[2]
  y.grid <- seq(min(y)*.66, max(y)*1.5, length = G)
  if (histogram) {
    hist(y, freq = FALSE, breaks = nbin, col = gray(.5), border = gray(.25), xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "")
  }
  else {
    plot(c(0, y.grid), c(0, apply(fdens, 2, mean)), type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  }
  if (bands) {
    color <- rgb(190, 190, 190, alpha=180, maxColorValue=255)
    grid <- c(y.grid, rev(y.grid))
    perc <- (1 - confid) / 2
    curves <- c(allcurves.q(fdens, perc), rev(allcurves.q(fdens, (1 - perc))))
    polygon(grid, curves, col = color, lty = 1, lwd = 2, border = NA)
  }
  for (i in sample(1:dim(fdens)[1], ndens, replace = rpl)) lines(y.grid, fdens[i, ], col = "red")
  lines(y.grid, apply(fdens, 2, mean), lwd = 2)
  rug(y)
}

