\name{prob.predict}
\alias{prob.predict}
\title{Tail probability estimation for a Gamma Shape Mixture Model}
\usage{
   prob.predict(mcmc.w,mcmc.theta,thresh)
}
\description{
   This function allows to estimate the tail probability of a Gamma Shape Mixture Model using the output of the \code{\link{gsm}} or \code{\link{gsm.theta}} procedures.
}
\arguments{
   \item{mcmc.w}{matrix of the mixture's weights posterior draws; it is part of the output of the \code{\link{gsm}} or \code{\link{gsm.theta}} functions.}
   \item{mcmc.theta}{vector of the mixture's rate parameter posterior draws; it is part of the output of the \code{\link{gsm}} or \code{\link{gsm.theta}} functions.}
   \item{thresh}{threshold value.}
}
\details{
   The tail probability is estimated by applying the standard Rao-Blackwellized estimator on the Gibbs sampling realizations obtained through the \code{\link{gsm}} or \code{\link{gsm.theta}} procedures.
}
\value{
   A numerical vector containing the posterior draws for the tail probability exceeding the value of \code{thresh}.
}
\references{
   Venturini, S., Dominici, F., and Parmigiani, G., "Gamma Shape Mixtures for Heavy-Tailed Distributions" (December 2006). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper 124. 
   \url{http://www.bepress.com/jhubiostat/paper124}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{gsm}},
   \code{\link{gsm.plot}}.
}
\keyword{distribution}
\keyword{models}
\examples{
set.seed(2040)
y <- rgsm(500,c(.1,.3,.4,.2),1)
burnin <- 100
J <- 250
gsm.out <- gsm(y,J,300,burnin+500,6500,340,1/J)
thresh <- c(0.1,0.5,0.75,1,2)
tail.prob.est <- rep(NA,length(thresh))
tail.prob.true <- rep(NA,length(thresh))
for (i in 1:length(thresh)){
   tail.prob.est[i] <- mean(prob.predict(gsm.out$weight[(burnin+1):600,],gsm.out$theta[(burnin+1):600],thresh[i]))
   tail.prob.true[i] <- sum(y>thresh[i])/length(y)
}
qqplot(tail.prob.true,tail.prob.est,main="Q-Q plot of true vs. estimated tail probability")
abline(0,1,lty=2)
}
