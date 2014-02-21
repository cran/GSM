\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,gsm,missing-method}
\title{Plot of a Gamma Shape Mixture Model}
\usage{
\S4method{plot}{gsm,missing}(x, ndens = 5, xlab = "x", ylab = "density", nbin = 10,
		histogram = FALSE, bands = FALSE, confid = .95, start = 1, \dots)
}
\description{
   \code{plot} method for class "gsm". This function plots the output of a Gamma Shape Mixture estimation procedure.
}
\arguments{
   \item{x}{object of class "gsm"; a list returned by the \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}} functions.}
   \item{ndens}{number of simulated density curves to plot.}
   \item{xlab}{a title for the x axis.}
   \item{ylab}{a title for the y axis.}
   \item{nbin}{number of bins for the histogram.}
   \item{histogram}{logical; if TRUE the histogram is plotted on the figure.}
   \item{bands}{logical; if TRUE the 95\% credibility bands are overimposed on the density graph.}
   \item{confid}{confidence level for the pointwise credibility bands around the density estimate.}
   \item{start}{MCMC run to start from.}
   \item{\dots}{further arguments passed to or from other methods.}
}
\details{
   To produce a standard histogram with the estimated density curve superimposed on it, simply set \code{ndens} to 0 and \code{histogram} to \code{TRUE}.
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }

\seealso{
   \code{\link{estim.gsm}},
   \code{\link{estim.gsm_theta}},
   \code{\link{summary-methods}},
   \code{\link{predict-methods}}.
}
\value{
   List with the following components:
  \item{xval}{horizontal coordinates.}
  \item{yval}{vertical coordinates (pointwise density posterior means).}
}
\keyword{methods}
\examples{
set.seed(2040)
y <- rgsm(500, c(.1, .3, .4, .2), 1)
burnin <- 5
mcmcsim <- 10
J <- 250
gsm.out <- estim.gsm(y, J, 300, burnin + mcmcsim, 6500, 340, 1/J)
par(mfrow = c(3, 2))
plot(gsm.out)
plot(gsm.out, ndens = 0, nbin = 20, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, histogram = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, histogram = TRUE, bands = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 5, nbin = 20, histogram = TRUE, bands = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, bands = TRUE, start = (burnin + 1))
}
