\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,gsm-method}
\title{Summarizing Gamma Shape Mixtures}
\description{
	\code{summary} method for class "gsm". This function allows to summarize the output of a Gamma Shape Mixture estimate procedure like \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}}.
}
\usage{
\S4method{summary}{gsm}(object, plot = FALSE, start = 1, \dots)
}
\arguments{
   \item{object}{object of class "gsm"; a list returned by the \code{\link{estim.gsm}} or \code{\link{estim.gsm_theta}} functions.}
   \item{plot}{logical; if TRUE produces a bar plot of the mixture weights posterior means.}
   \item{start}{MCMC run to start from.}
   \item{\dots}{further arguments passed to or from other methods.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
	The function \code{summary} computes and returns a list of summary statistics of the fitted gamma shape mixture given in \code{object}, in particular
   \item{theta}{summary index of the theta parameter posterior draws.}
   \item{weight posterior means}{vector of the mixture weights posterior means.}
}
\references{
   Venturini, S., Dominici, F. and Parmigiani, G. (2008), "Gamma shape mixtures for heavy-tailed distributions". Annals of Applied Statistics, \bold{Volume 2}, Number 2, 756--776.
   \url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{estim.gsm}},
   \code{\link{estim.gsm_theta}},
   \code{\link{plot-methods}},
   \code{\link{predict-methods}}.
}
\keyword{methods}
\examples{
set.seed(2040)
y <- rgsm(500, c(.1, .3, .4, .2), 1)
burnin <- 5
mcmcsim <- 10
J <- 250
gsm.out <- estim.gsm(y, J, 300, burnin + mcmcsim, 6500, 340, 1/J)
summary(gsm.out, TRUE, start = (burnin + 1))
}
