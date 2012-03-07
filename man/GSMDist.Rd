\name{GSMDist}
\alias{GSMDist}
\alias{dgsm}
\alias{pgsm}
\alias{rgsm}
\alias{qgsm}
\title{Utility Function}
\usage{
   dgsm(x, weight, rateparam)
   pgsm(q, weight, rateparam, lower.t = TRUE)
   rgsm(n, weight, rateparam)
   qgsm(p, x = NULL, weight, rateparam, alpha = .05, br = c(0, 1000), lower.t = TRUE)
}
\description{
   Function evaluations for a Gamma Shape Mixture Model.
}
\arguments{
   \item{x, q}{vector of quantiles.}
   \item{n}{number of observations.}
   \item{p}{vector of probabilities.}
   \item{weight}{vector of mixture weights.}
   \item{rateparam}{reciprocal of the shape parameter, as in \code{\link{GammaDist}}.}
   \item{alpha}{outside the interval (alpha, 1 - alpha) the quantiles are found by searching for the root of \emph{F(x) - p = 0}.}
   \item{br}{a vector containing the end-points of the interval to be searched for the root.}
   \item{lower.t}{logical; if TRUE (default), probabilities are \emph{P[X <= x]} otherwise, \emph{P[X > x]}.}
}
\details{
   The parametrisation implemented in this function is described in Venturini et al. (2008).
}
\value{
   \code{\link{dgsm}} gives the density, \code{\link{pgsm}} gives the distribution function, \code{\link{qgsm}} gives the quantile function, and \code{\link{rgsm}} generates random deviates.
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\seealso{
   \code{\link{dgamma}},
   \code{\link{pgamma}},
   \code{\link{rgamma}},
   \code{\link{uniroot}}.
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\keyword{distribution}
