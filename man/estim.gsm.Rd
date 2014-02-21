\name{estim.gsm}
\alias{estim.gsm}
\title{Estimation of a Gamma Shape Mixture Model (GSM) with collapsing}
\description{
   This function provides the inferential algorithm to estimate a mixture of gamma distributions in which the mixing occurs over the shape parameter. It implements the collapsing approach for the GSM model, as discussed in Venturini et al. (2008).
}
\usage{
   estim.gsm(y, J, G = 100, M = 600, a, b, alpha, init = list(rep(1 / J, J), NA,
   		rep(1, N)))
}
\arguments{
   \item{y}{vector of data.}
   \item{J}{number of mixture components.}
   \item{G}{number of points where to evaluate the GSM density.}
   \item{M}{number of MCMC runs.}
   \item{a}{hyperparameter of the rate parameter prior distribution.}
   \item{b}{hyperparameter of the rate parameter prior distribution.}
   \item{alpha}{hyperparameter of the mixture's weights prior distribution.}
   \item{init}{initialization values.}
}
\details{
   Suggestions on how to choose \code{J}, \code{a} and \code{b} are provided in Venturini et al. (2008). In that work the \code{alpha} vector is always set at (1/\code{J},...,1/\code{J}), but here one is free to choose the value of the generic element of \code{alpha}.
}
\value{
   \code{estim.gsm} returns an object of \link{class} "gsm", which is a list with the following components:
   \item{fdens}{matrix containing the posterior draws for the mixture's density.}
   \item{theta}{vector containing the posterior draws for the mixture's rate parameter.}
   \item{weight}{matrix containing the posterior draws for the mixture's weights.}
   \item{label}{matrix containing the posterior draws for the mixture's labels.}
   \item{data}{vector of data.}
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{estim.gsm_theta}}, 
   \code{\link{summary-methods}},
   \code{\link{plot-methods}}.
}
\keyword{distribution}
\keyword{mixture}
\examples{
\dontrun{
set.seed(2040)
y <- rgsm(500, c(.1, .3, .4, .2), 1)
burnin <- 100
mcmcsim <- 500
J <- 250
gsm.out <- estim.gsm(y, J, 300, burnin + mcmcsim, 6500, 340, 1/J)
summary(gsm.out, plot = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, histogram = TRUE, start = (burnin + 1))}
}
