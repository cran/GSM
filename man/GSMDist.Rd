\name{GSMDist}
\alias{GSMDist}
\alias{dgsm}
\alias{pgsm}
\alias{rgsm}
\title{Utility Function}
\usage{
   dgsm(x, weight, rateparam)
   pgsm(q, weight, rateparam)
   rgsm(n, weight, rateparam)
}
\description{
   Density function evaluation for a Gamma Shape Mixture Model.
}
\arguments{
   \item{x, q}{vector of quantiles.}
   \item{n}{number of observations.}
   \item{weight}{vector of mixture weights.}
   \item{rateparam}{reciprocal of the shape parameter, as in \code{\link{GammaDist}}.}
}
\details{
   The parametrisation implemented in this function is described in Venturini et al. (2006).
}
\value{
   \code{\link{dgsm}} gives the density, \code{\link{pgsm}} gives the distribution function, and \code{\link{rgsm}} generates random deviates.
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\seealso{
   \code{\link{dgamma}},
   \code{\link{pgamma}},
   \code{\link{rgamma}}.
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\keyword{distribution}
