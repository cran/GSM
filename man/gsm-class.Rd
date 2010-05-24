\name{gsm-class}
\Rdversion{1.1}
\docType{class}
\alias{gsm-class}
\alias{initialize,gsm-method}
%\alias{summary,gsm-method}
%\alias{plot,gsm-method}
%\alias{predict,gsm-method}

\title{Class "gsm". Result of Gamma Shape Mxiture Estimation.}
\description{
  This class encapsulates results of a Gamma Shape Mixture estimation procedure.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gsm", fdens, theta, weight, data)}, but most often as the result of a call to \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}}.
}
\section{Slots}{
  \describe{
    \item{\code{fdens}:}{Object of class \code{"matrix"}; posterior draws from the MCMC simulation algorithm of the Gamma Shape Mixture density.}
    \item{\code{theta}:}{Object of class \code{"numeric"}; posterior draws from the MCMC simulation algorithm of the Gamma Shape Mixture scale parameter.}
    \item{\code{weight}:}{Object of class \code{"matrix"}; posterior draws from the MCMC simulation algorithm of the Gamma Shape Mixture weights.}
    \item{\code{label}:}{Object of class \code{"matrix"}; posterior draws from the MCMC simulation algorithm of the Gamma Shape Mixture lables.}
    \item{\code{data}:}{Object of class \code{"numeric"}; original data.}
  }
}
%\section{Extends}{
%Class \code{"\linkS4class{list}"}, from data part.
%Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
%}
\section{Methods}{
  \describe{
%    \item{initialize}{\code{signature(.Object = "gsm")}: Class initialization.}
    \item{plot}{\code{signature(x = "gsm", y = "missing")}: Plot Gamma Shape Mixture estimate.}
    \item{predict}{\code{signature(object = "gsm")}: Estimate of the Gamma Shape Mixture upper tail.}
    \item{summary}{\code{signature(object = "gsm")}: Generate object summary.}
	 }
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{estim.gsm}}, 
   \code{\link{summary-methods}},
   \code{\link{plot-methods}},
   \code{\link{predict-methods}},
   \code{\link{summary-methods}}.
}
\keyword{methods}
\keyword{class}
