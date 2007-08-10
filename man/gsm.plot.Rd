\name{gsm.plot}
\alias{gsm.plot}
\title{Plot of a Gamma Shape Mixture Model}
\description{
   This function plots the output of a Gamma Shape Mixture estimation procedure.
}
\usage{
   gsm.plot(v,y,ndens=5,xlim=c(min(y),max(y)),ylim=c(0,max(v$ff)),xlab="x",ylab="density",nbin=10,histogram=FALSE,bands=FALSE)
}
\arguments{
   \item{v}{list returned by the \code{\link{gsm}} or \code{\link{gsm.theta}} functions.}
   \item{y}{vector of data to use in the analysis.}
   \item{ndens}{number of simulated density curves to plot.}
   \item{xlim}{the x limits (x1, x2) of the plot.}
   \item{ylim}{the y limits of the plot.}
   \item{xlab}{a title for the x axis.}
   \item{ylab}{a title for the y axis.}
   \item{nbin}{number of bins for the histogram.}
   \item{histogram}{logical; if TRUE the histogram is plotted on the figure.}
   \item{bands}{logical; if TRUE the 95\% credibility bands are overimposed on the density graph.}
}
\details{
   To produce a standard histogram with the estimated density curve superimposed on it, simply set \code{ndens} to 0 and \code{histogram} to \code{TRUE}.
}
\references{
   Venturini, S., Dominici, F., and Parmigiani, G., "Gamma Shape Mixtures for Heavy-Tailed Distributions" (December 2006). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper 124. 
   \url{http://www.bepress.com/jhubiostat/paper124}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }

\seealso{
   \code{\link{gsm}}, 
   \code{\link{gsm.theta}}, 
   \code{\link{rgsm}}.
}
\keyword{distribution}
\keyword{models}
\examples{
set.seed(2040)
y <- rgsm(500,c(.1,.3,.4,.2),1)
burnin <- 100
gsm.out <- gsm(y,250,300,burnin+500,6500,340,.25)
par(mfrow=c(2,2))
gsm.plot(gsm.out,y,ndens=5,nbin=20)
gsm.plot(gsm.out,y,ndens=0,nbin=20,histogram=TRUE)
gsm.plot(gsm.out,y,ndens=0,nbin=20,histogram=TRUE,bands=TRUE)
gsm.plot(gsm.out,y,ndens=5,nbin=20,histogram=TRUE,bands=TRUE)
}
