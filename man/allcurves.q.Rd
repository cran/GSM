\name{allcurves.q}
\alias{allcurves.q}
\title{Utility function}
\usage{
   allcurves.q(postdata,perc)
}
\description{
   Utility function for plotting a Gamma Shape Mixture Model density.
}
\arguments{
   \item{postdata}{matrix containing of a mixture's density posterior draws.}
   \item{perc}{percentile.}
}
\details{This is a utility function used to generate the credibility bands for a Gamma Shape Mixture density within \code{\link{gsm.plot}}.}
\seealso{
   \code{\link{gsm.plot}}.
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\keyword{distribution}
