\name{gsm}
\alias{gsm}
\title{Estimation of a Gamma Shape Mixture Model (GSM)}
\description{
   This function provides the inferential algorithm to estimate a mixture of gamma distributions in which the mixing occurs over the shape parameter. It implements the collapsing approach for the GSM model, as discussed in Venturini et al. (2006).
}
\usage{
   gsm(y,J,G,M,a,b,alpha)
}
\arguments{
   \item{y}{vector of data.}
   \item{J}{number of mixture components.}
   \item{G}{number of points where to evaluate the GSM density.}
   \item{M}{number of MCMC runs.}
   \item{a}{hyperparameter of the rate parameter prior distribution.}
   \item{b}{hyperparameter of the rate parameter prior distribution.}
   \item{alpha}{hyperparameter of the mixture's weights prior distribution.}
}
\details{
   Suggestions on how to choose \code{J}, \code{a} and \code{b} are provided in Venturini et al. (2006). In that work the \code{alpha} vector is always set at (1/\code{J},...,1/\code{J}), but here one is free to choose the value of the generic element of \code{alpha}.
}
\value{
   List with the following components:
   \item{J}{number of mixture components used in the GSM model.}
   \item{a}{hyperparameter of the rate parameter prior distribution used in the GSM model.}
   \item{b}{hyperparameter of the rate parameter prior distribution used in the GSM model.}
   \item{alpha}{hyperparameter of the mixture's weights prior distribution used in the GSM model.}
   \item{ff}{matrix containing the posterior draws for the mixture's density.}
   \item{y.grid}{vector of values used to evalute the GSM density.}
   \item{theta}{vector containing the posterior draws for the mixture's rate parameter.}
   \item{label}{matrix containing the posterior draws for the mixture's hidden label.}
   \item{weight}{matrix containing the posterior draws for the mixture's weights.}
}
\references{
   Venturini, S., Dominici, F., and Parmigiani, G., "Gamma Shape Mixtures for Heavy-Tailed Distributions" (December 2006). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper 124. 
   \url{http://www.bepress.com/jhubiostat/paper124}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{gsm.theta}}, 
   \code{\link{gsm.plot}}.
}
\keyword{distribution}
\keyword{models}
\examples{
set.seed(2040)
y <- rgsm(500,c(.1,.3,.4,.2),1)
burnin <- 100
gsm.out <- gsm(y,250,300,burnin+500,6500,340,.25)
gsm.plot(gsm.out,y,ndens=0,nbin=20,histogram=TRUE)
}
