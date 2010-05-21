\name{allcurves.q}
\alias{allcurves.q}
\title{Utility Function}
\usage{
   allcurves.q(post, perc)
}
\description{
   Utility function for plotting a Gamma Shape Mixture Model density.
}
\arguments{
   \item{post}{matrix containing of a mixture's density posterior draws.}
   \item{perc}{percentile, a value that satisfies 0 < perc < 1.}
}
\details{This is a utility function used to generate the credibility bands for a Gamma Shape Mixture density within \code{\link{plot}}.}
\seealso{
   \code{\link{plot-methods}}.
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
%\keyword{distribution}
