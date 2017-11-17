#' @title Rank abundance distribution curve
#' @description Given a community abundance matrix, compute the species abundance distribution as rank abundance curve.
#' @details Please provide species abundances as counts (scale and round if necessary). Note that for fit.rad, Null is the broken stick distribution.
#' A vector of species abundances can be provided as well.
#' @param x community matrix with species as rows, samples as columns and non-negative integers as entries
#' @param remove.zeros if true, zeros in the species abundance distribution curve are discarded
#' @param fit.neutral assess the fit of the neutral model with function theta.prob from the untb R package
#' @param fit.rad (optional) fit different RAD models using vegan's radfit function. If selected, fit.distrib is ignored and if plot is true, vegan's plot function for RAD curves is used.
#' @param fit.distrib (optional) select the best-fitting distribution for the RAD by comparing normal, lognormal, gamma, negative binomial and Poisson distribution using MASS function fitdistr
#' @param fit.N (optional) number of samples drawn from the distribution to be fitted (defaults to the number of species in x after filtering of zeros)
#' @param sample.index (optional) the sample index for which species abundance distribution is computed, if not provided, the mean across samples is taken
#' @param plot plot the RAD curve
#' @param header add the header to the title of the plot
#' @return the RAD parameters are returned, if fit.neutral is true, the probability of the biodiversity parameter theta is returned, if fit.rad or fit.distrib is true, the results of distribution fitting are also returned (distrib, score [log-likelihood for fit.distrib and AIC for fit.rad], params and fitted.rad)
#' @examples
#' ts.ricker=ricker(10,generateA(10),K=rep(0.01,10))
#' ts.ricker=round(ts.ricker*1000) # scale to integers
#' rad.out.ricker=rad.fit.ricker=rad(ts.ricker, fit.distrib = TRUE, header="Ricker model")
#' ts.hubbell=simHubbell(N = 100, M = 100, I = 500, d = 1, tskip=500, tend=1500)
#' rad.out.hubbell=rad.fit.hubbell=rad(ts.hubbell, fit.distrib = TRUE, header="neutral model")
#' @export

rad<-function(x,remove.zeros=TRUE, fit.neutral=FALSE, fit.rad=FALSE, fit.distrib=FALSE,fit.N=nrow(x), sample.index=NA, plot=TRUE, header=""){
  # identify inflection points: library(bcp)
  # chps=bcp(sad.out$sad)
  # which(chps$posterior.prob==max(chps$posterior.prob,na.rm=TRUE))
  if(!is.matrix(x)){
    x=matrix(x,nrow=length(x),ncol=1)
  }
  distrib=""
  score=NA
  thetaprob=NA
  params=c()
  fitted.sad=c()
  rad.out=NULL
  if(fit.neutral==TRUE && (fit.rad==TRUE || fit.distrib==TRUE)){
    stop("Fitting neutral model excludes fit.rad and fit.distrib.")
  }
  if(!is.whole.matrix(x)){
    stop("The matrix should consist of integers! Please scale.")
  }
  if(!is.na(sample.index)){
    vec=x[,sample.index]
  }else if(ncol(x)>1){
    vec=round(apply(x,1,mean))
  }else{
    vec=x[,1]
  }
  if(fit.rad==TRUE && fit.distrib==TRUE){
    warning("Can either do fit.rad or fit.distrib. Continuing with fit.rad.")
    fit.distrib==FALSE
  }
  if(remove.zeros==TRUE){
    zero.indices=which(vec==0)
    if(length(zero.indices)>0){
      print(paste("Removing",length(zero.indices),"species with zero abundance."))
      ok.indices=setdiff(1:length(vec),zero.indices)
      vec=vec[ok.indices]
      if(fit.N==nrow(x)){
        fit.N=length(vec)
      }
    }
  }
  vec=sort(vec,decreasing = TRUE)
  if(fit.distrib==TRUE){
    fit.res=identifyDistrib(vec, N=fit.N)
    distrib=fit.res$distrib
    params=fit.res$params
    score=fit.res$loglik
    fitted.sad=fit.res$sad
    #print(fit.res$sad)
  }else if(fit.rad==TRUE){
    rad.out=radfit(vec)
    distribs=c("Null","Zipf","Preemption","Mandelbrot","Lognormal")
    akaikes=c(rad.out$models$Null$aic,rad.out$models$Zipf$aic,rad.out$models$Preemption$aic,rad.out$models$Mandelbrot$aic,rad.out$models$Lognormal$aic)
    # if some distributions got the same log likelihood, choose one at random
    akaikes=akaikes+runif(length(akaikes),min=0,max=0.001)
    winner.index=which(akaikes==min(akaikes))
    distrib=distribs[winner.index]
    score=akaikes[winner.index]
    params=rad.out$models[[distrib]]$coefficients
    fitted.sad=rad.out$models[[distrib]]$fitted.values
  }else if(fit.neutral==TRUE){
    thetaprob=untb::theta.prob(theta=untb::optimal.theta(vec),x=vec,give.log=FALSE)
  }
  if(fit.rad==TRUE || fit.distrib==TRUE){
    print(paste("Selected distribution: ",distrib))
    print(paste("Fitting score: ",score))
    print("Fitting parameters:")
    print(params)
  }
  # TODO: generate neutral data with untb and plot together with RAD
  if(plot==TRUE){
    if(fit.rad==TRUE){
      plot(rad.out)
    }else{
      col1=rgb(0,1,0,0.5)
      col2=rgb(1,0,0,0.5)
      ylim=c(0,max(vec))
      if(fit.distrib==TRUE){
        ylim=c(0,max(max(vec),max(fitted.sad)))
        if(length(fitted.sad) > length(vec)){
          # add dummy bars
          diff=length(fitted.sad)-length(vec)
          vec=c(vec,rep(NA,diff))
        }
      }
      barplot(vec, xlab="Rank", ylab="Abundance", main=paste("RAD",header), ylim=ylim, col=col1)
      #abline()
      if(fit.distrib==TRUE){
        barplot(fitted.sad,col=col2, add=TRUE)
        legend("topright",legend=c("observed",fit.res$distrib), lty = rep(1,2), col = c(col1,col2), merge = TRUE, bg = "white", text.col="black")
      }
    }
  }
  res=list()
  res[["rad"]]=vec
  if(fit.distrib==TRUE || fit.rad==TRUE){
    res[["distrib"]]=distrib
    res[["score"]]=score
    res[["params"]]=params
    res[["fitted.rad"]]=fitted.sad
  }else if(fit.neutral==TRUE){
    res[["thetaprob"]]=thetaprob
  }
  return(res)
}


# Identify the distribution for a given vector using MASS::fitdistr
# Tested distributions include normal, lognormal, gamma, negative binomial and Poisson
# If generateSAD is true, a species-abundance vector with N entries is produced with the best-fitting distribution.
identifyDistrib<-function(v, N=length(v), generateSAD=TRUE){
  distribtypes=c("normal","lognormal","negative binomial", "gamma", "Poisson")
  win.loglik=NA
  win.distrib="NA"
  win.params=list()
  failure=FALSE
  sad=c()
  if(sd(v)>0){
    logliks=c()
    params=list()
    index=1
    # loop supported distributions
    for(distribname in distribtypes){
      fit.out=NULL
      # following: http://mazamascience.com/WorkingWithData/?p=912
      res<-tryCatch({
        fit.out=MASS::fitdistr(v,densfun=distribname)
      },
      warning=function(war){print(paste("WARN:",war)); return(NULL)},
      error=function(err){print(paste("ERROR:",err)); return(NULL)}
      )
      if(!is.null(fit.out)){
        logliks=c(logliks,fit.out$loglik)
        params[[index]]=fit.out$estimate
        #print(distribname)
        #print(fit.out$estimate)
      }else{
        logliks=c(logliks,NA)
        params[[index]]=NA
        failure=TRUE
        print(paste("Fitting of distribution",distribname,"failed."))
      }
      index=index+1
    } # loop distributions
    # if some distributions got the same log likelihood, choose one at random
    logliks=logliks+runif(length(logliks),min=0,max=0.001)
    max.index=which(logliks==max(logliks,na.rm=TRUE))
    if(length(max.index) > 0){
      win.distrib=distribtypes[max.index]
      win.loglik=logliks[max.index]
      win.params=params[[max.index]]
    }
    # generate SAD curve with winning distribution
    if(generateSAD==TRUE){
      if(win.distrib=="normal"){
        sad=sort(rnorm(N,mean=win.params[1],sd=win.params[2]), decreasing=TRUE)
      }else if(win.distrib=="Poisson"){
        sad=sort(rpois(N,lambda=win.params[1]), decreasing=TRUE)
      }else if(win.distrib=="gamma"){
        sad=sort(rgamma(N,shape=win.params[1],rate=win.params[2]),decreasing=TRUE)
      }else if(win.distrib=="lognormal"){
        sad=sort(rlnorm(N,meanlog=win.params[1],sdlog=win.params[2]),decreasing=TRUE)
      }else if(win.distrib=="negative binomial"){
        sad=sort(rnbinom(N,size=win.params[1], mu=win.params[2]),decreasing=TRUE)
      }
    }
  } # v is not all zero
  else{
    warning("Vector is constant!")
  }

  res=list(win.distrib, win.loglik, win.params, sad, failure)
  names(res)=c("distrib","loglik", "params", "sad","failurehappened")
  return(res)
}

# Identify distributions in a matrix row-wise.
#
# param x a count matrix
# param pseudocount how to treat zeros: if false, omit species with zero abundance, if true, add a pseudocount of one
# param skipFailures: do not assign a distribution when fitting one of the test distributions failed
#
# examples
# tsR=ricker(N=100,A=generateA(N=100),y=generateAbundances(N=100, mode=5, probabs=TRUE),tend=500)
# countsR=round(ts[,100:ncol(tsR)]*1000) # skip transient and convert to counts
# distribsR=identifyDistribs(countsR)
#
# tsN=rbind(rnorm(100,mean=10),rnorm(100,mean=10))
# countsN=round(tsN*1000) # convert to counts
# distribsN=identifyDistribs(countsN)
#
# tsH=simHubbell(N=10,M=100,I=1500,d=10,tend=1000,tskip=500,m.vector=generateAbundances(N=100, mode=5, probabs=TRUE))
# distribsH=identifyDistribs(tsH)
#
# takes ages:
# tsS=soi(100, I=1500, A=generateA(N=100,c=0.01), m.vector=generateAbundances(N=100,mode=5), tend=100)
# distribsS=identifyDistribs(tsS[,20:ncol(tsS)])

identifyDistribs<-function(x, pseudocount=FALSE, skipFailures=TRUE){
  distribs=c()
  for(i in 1:nrow(x)){
    v=as.numeric(x[i,])
    zero.indices=which(v==0)
    if(length(zero.indices)>0){
      # remove zeros
      if(pseudocount==FALSE){
        ok.indices=setdiff(1:length(v),zero.indices)
        v=v[ok.indices]
      }else{
        v[zero.indices]=1
      }
    }
    if(sd(v)>0){
     fit.res=identifyDistrib(v, generateSAD = FALSE)
     if(fit.res$failurehappened==TRUE && skipFailures==TRUE){
       distribs=c(distribs,"NA")
     }else{
      distribs=c(distribs,fit.res$distrib)
     }
    }else{
      distribs=c(distribs,"NA")
    }
  } # end loop rows
  distribsummary=list(which(distribs=="normal"),which(distribs=="lognormal"),which(distribs=="negative binomial"),which(distribs=="Poisson"),which(distribs=="gamma"),which(distribs=="NA"))
  names(distribsummary)=c("norm","lognorm","negbin","pois","gamma","none")
  print(paste("Number of taxa with normal distribution: ",length(distribsummary$norm)))
  print(paste("Number of taxa with lognormal distribution: ",length(distribsummary$lognorm)))
  print(paste("Number of taxa with negative binomial distribution: ",length(distribsummary$negbin)))
  print(paste("Number of taxa with Poisson distribution: ",length(distribsummary$pois)))
  print(paste("Number of taxa with Gamma distribution: ",length(distribsummary$gamma)))
  print(paste("Not successfully fitted to tested distributions: ",length(distribsummary$none)))
  return(distribsummary)
}


# Plot histograms of observed and estimated abundances.
# Abundances are estimated with a random distribution.
#
# Note that density estimation with function density does not work
# well. We also need to manually tune the adjust value
# (6 was OK for test time series) and the band width (default
# did not work well, but SJ is too slow for 100x3000 matrices).
#
# x: matrix with species as rows and time points as columns
# distrib: distribution used to estimate abundances, can be pois (Poisson), norm (Gaussian) or negbin (negative binomial)
# data: a string describing the data type for the plot title
#
# Example:
# fitDistrib(timeseries$exp36, data="SOI 36")
fitDistrib<-function(x, distrib="pois", data=""){
  pooled=as.numeric(x)
  mean=mean(pooled)
  var=var(pooled)
  N=length(pooled)
  distribName=""
  print(paste("mean:",mean))
  print(paste("var:",var))
  print(paste("Distribution:",distrib))
  if(distrib=="pois"){
    sim=rpois(N,mean)
    distribName="Poisson"
  }else if(distrib=="norm"){
    sim=rnorm(N,mean=mean,sd=sqrt(var))
    distribName="Gaussian"
  }else if(distrib=="negbin"){
    sim=rpois(N,rgamma(N,mean))
    distribName="negative Binomial"
  }
  col1=rgb(0,1,0,0.5)
  col2=rgb(1,0,0,0.5)
  # limits
  xmax=max(pooled,na.rm=TRUE)
  xmin=min(pooled,na.rm=TRUE)
  ymax=max(sim,na.rm=TRUE)
  ymin=min(sim,na.rm=TRUE)
  max=max(xmax,ymax)
  min=min(ymin,xmin)
  xmaxD=max(hist(pooled,breaks="FD",plot=FALSE)$density)
  ymaxD=max(hist(sim,breaks="FD",plot=FALSE)$density)
  maxD=max(xmaxD,ymaxD)
  maxD=maxD+0.5 # add a margin
  title=paste("Histogram of observed abundances and \nabundances estimated with ",distribName," for ",data,sep="")
  hist(pooled,breaks="FD",xlim=c(min,max), ylim=c(0,maxD), prob=TRUE,col=col1, border=col1,xlab="Abundances", main=title)
  hist(sim,breaks="FD",prob=TRUE,col=col2, border=col2,add=TRUE)
  legend("right",legend=c("Observed","Estimated"), lty = rep(1,2), col = c(col1,col2), merge = TRUE, bg = "white", text.col="black")
}
