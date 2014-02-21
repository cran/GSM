\name{predict-methods}
\docType{methods}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,gsm-method}
\title{Tail Probability Estimation for a Gamma Shape Mixture Model}
\usage{
\S4method{predict}{gsm}(object, thresh, start = 1, \dots)
}
\description{
   \code{predict} method for class "gsm". This function allows to estimate the tail probability of a Gamma Shape Mixture Model using the output of the \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}} procedures.
}
\arguments{
   \item{object}{object of class "gsm"; a list returned by the \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}} functions.}
   \item{thresh}{threshold value.}
   \item{start}{MCMC run to start from.}
   \item{\dots}{further arguments passed to or from other methods.}
}
\details{
   The tail probability is estimated by applying the standard Rao-Blackwellized estimator on the Gibbs sampling realizations obtained through the \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}} procedures.
}
\value{
   A numerical vector containing the posterior draws for the tail probability exceeding the value of \code{thresh}.
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{estim.gsm}},
   \code{\link{estim.gsm_theta}},
   \code{\link{predict-methods}},
   \code{\link{plot-methods}}.
}
\keyword{methods}
\examples{
set.seed(2040)
y <- rgsm(500, c(.1, .3, .4, .2), 1)
burnin <- 5
mcmcsim <- 10
J <- 250
gsm.out <- estim.gsm(y, J, 300, burnin + mcmcsim, 6500, 340, 1/J)
thresh <- c(0.1, 0.5, 0.75, 1, 2)
tail.prob.est <- tail.prob.true <- rep(NA, length(thresh))
for (i in 1:length(thresh)){
   tail.prob.est[i] <- mean(predict(gsm.out, thresh[i]))
   tail.prob.true[i] <- sum(y > thresh[i])/length(y)
}
qqplot(tail.prob.true, tail.prob.est, main = "Q-Q plot of true vs. estimated tail probability")
abline(0, 1, lty = 2)
}
