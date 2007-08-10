\name{GSMDist}
\alias{GSMDist}
\alias{dgsm}
\alias{pgsm}
\alias{rgsm}
\title{Utility function}
\usage{
   dgsm(x,weight,rateparam)
   pgsm(q,weight,rateparam)
   rgsm(n,weight,rateparam)
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
   Venturini, S., Dominici, F., and Parmigiani, G., "Gamma Shape Mixtures for Heavy-Tailed Distributions" (December 2006). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper 124. 
   \url{http://www.bepress.com/jhubiostat/paper124}
}
\seealso{
   \code{\link{dgamma}},
   \code{\link{pgamma}},
   \code{\link{rgamma}}.
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\keyword{distribution}
