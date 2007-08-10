gsm <- function(y,J,G,M,a,b,alpha){
  N <- length(y)
  y.grid <- seq(min(y)*.66,max(y)*1.5,length=G)

  lbl <- matrix(NA,nrow=N,ncol=M)
  wgt <- matrix(NA,nrow=M,ncol=J)

  x <- rep(1,N)
  w <- rep(1/J,J)
  theta <- rep(NA,M)
  logff <- matrix(NA,J,G)
  ff <- matrix(NA,M,G)

  for (m in 1:M){
    for (nn in 1:N) {
      temp <- log(seq(a+sum(x)-x[nn],a+sum(x)-x[nn]+J-1))
      pi <- (1:J-1)*log(y[nn]) - lgamma(1:J) + cumsum(temp) - (1:J)*log(b+sum(y))
      pi <- ( w*exp(pi) ) / sum( w*exp(pi) )
      x[nn] <- sample(1:J,1,prob=pi)
      lbl[nn,m]=x[nn]
    }

    x.counts <- table( factor(x, levels=1:J ))
    w <- rdirichlet(1, x.counts + alpha )
    theta[m] <- rgamma( 1, sum(x)+a, rate=sum(y)+b )
    for (j in 1:J){
      logff[j,] <- log(w[j]) + (j-1)*log(y.grid) - theta[m]*y.grid - lgamma(j) + j*log(theta[m])
    }
    ff[m,] <- apply( exp( logff ), 2, sum )
    wgt[m,] <- w
    if (m/100==round(m/100)) print(paste("simulation",m,"/",M))
  }
  return(list(J=J,a=a,b=b,alpha=alpha,ff=ff,y.grid=y.grid,theta=theta,label=lbl,weight=wgt))
}

gsm.theta <- function(y,J,G,M,a,b,alpha){
  N <- length(y)
  y.grid <- seq(min(y)*.66,max(y)*1.5,length=G)

  lbl <- matrix(NA,nrow=N,ncol=M)
  wgt <- matrix(NA,nrow=M,ncol=J)

  x <- rep(1,N)
  w <- rep(1/J,J)
  theta <- rep(J/max(y),M+1)
  logff <- matrix(NA,J,G)
  ff <- matrix(NA,M,G)

  for (m in 1:M){
    for (nn in 1:N) {
      pi <- (1:J-1)*log(y[nn]) - theta[m]*y[nn] - lgamma(1:J) + (1:J)*log(theta[m])
      pi <- ( w*exp(pi) ) / sum( w*exp(pi) )
      x[nn] <- sample(1:J,1,prob=pi)
      lbl[nn,m]=x[nn]
    }

    x.counts <- table( factor(x, levels=1:J ))
    w <- rdirichlet(1, x.counts + alpha )
    theta[m+1] <- rgamma( 1, sum(x)+a, rate=sum(y)+b )
    for (j in 1:J){
      logff[j,] <- log(w[j]) + (j-1)*log(y.grid) - theta[m+1]*y.grid -
         lgamma(j) + j*log(theta[m+1])
    }
    ff[m,] <- apply( exp( logff ), 2, sum )
    wgt[m,] <- w
    if (m/100==round(m/100)) print(paste("simulation",m,"/",M))
  }
  return(list(J=J,a=a,b=b,alpha=alpha,ff=ff,y.grid=y.grid,theta=theta,label=lbl,weight=wgt))
}

gsm.plot <- function(v,y,ndens=5,xlim=c(min(y),max(y)),ylim=c(0,max(v$ff)),xlab="x",ylab="density",nbin=10,histogram=FALSE,bands=FALSE){
  y.temp <- y
  if (histogram) {
    hist(y.temp,freq=FALSE,breaks=nbin,col="lightgray",border="gray",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main="")
    if (bands) {
        ygrid <- c(v$y.grid,rev(v$y.grid))
        bb <- c(allcurves.q(v$ff,0.025),rev(allcurves.q(v$ff,0.975)))
        polygon(ygrid,bb,col=gray(0.6),lty=1,lwd=2,border=NA)
    }
    lines(c(0,v$y.grid),c(0,apply(v$ff,2,mean)),type="n")
  }
  else {
    plot(c(0,v$y.grid),c(0,apply(v$ff,2,mean)),type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
    if (bands) {
        ygrid <- c(v$y.grid,rev(v$y.grid))
        bb <- c(allcurves.q(v$ff,0.025),rev(allcurves.q(v$ff,0.975)))
        polygon(ygrid,bb,col=gray(0.6),lty=1,lwd=2,border=NA)
    }
  }
  rug(y.temp)
  for(i in sample(1:dim(v$ff)[1],ndens)) lines(v$y.grid,v$ff[i,],col="red")
  lines(v$y.grid,apply(v$ff,2,mean),lwd=2)
}

allcurves.q <- function(postdata,perc){
  n <- dim(postdata)[2]
  temp <- rep(NA,n)
  for (i in 1:n) temp[i] <- quantile(postdata[,i],perc)
  return(temp)
}

prob.predict <- function(mcmc.w,mcmc.theta,thresh){
  numsim <- dim(mcmc.w)[1]
  pred.prob <- vector(length=numsim)
  for (i in 1:numsim) pred.prob[i] <- 1-pgsm(thresh,mcmc.w[i,],mcmc.theta[i])
  return(pred.prob)
}

dgsm <- function(x,weight,rateparam){
  numcomp <- length(weight)
  mixcomp <- matrix(NA,nrow=length(x),ncol=numcomp)
  for (i in 1:numcomp) mixcomp[,i] <- dgamma(x,shape=i,rate=rateparam)
  dens <- mixcomp%*%weight
  return(dens)
}

pgsm <- function(q,weight,rateparam){
  numcomp <- length(weight)
  mixcomp <- matrix(NA,nrow=length(q),ncol=numcomp)
  for (i in 1:numcomp) mixcomp[,i] <- pgamma(q,shape=i,rate=rateparam)
  cumprob <- mixcomp%*%weight
  return(cumprob)
}

rgsm <- function(n,weight,rateparam){
  J <- length(weight)
  rmixg <- vector(length=n)
  tmp.rmixg <- matrix(NA,nrow=J,ncol=n)
  tmp.lbl <- matrix(NA,nrow=J,ncol=n)
  for (i in 1:J) tmp.rmixg[i,] <- rgamma(n,shape=i,rate=rateparam)
  tmp.lbl <- rmultinom(n,1,weight)
  for (i in 1:J) rmixg <- rmixg + tmp.lbl[i,]*tmp.rmixg[i,]
  return(rmixg)
}
