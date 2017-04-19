# Identify distributions in a matrix row-wise.
#
# param x a count matrix
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
# tsS=soc(100, I=1500, A=generateA(N=100,c=0.01), m.vector=generateAbundances(N=100,mode=5), tend=100)
# distribsS=identifyDistribs(tsS[,20:ncol(tsS)])

identifyDistribs<-function(x){
  distribtypes=c("normal","lognormal","negative binomial", "Poisson")
  distribs=c()
  for(i in 1:nrow(x)){
    v=as.numeric(x[i,])
    if(sum(v)>0){
      logliks=c()
      # loop supported distributions
      for(distribname in distribtypes){
        #print(paste("Fitting",distribname,"to",i))
        # following: http://mazamascience.com/WorkingWithData/?p=912
        res<-tryCatch({
          fit.out=MASS::fitdistr(v,densfun=distribname)
        },
        warning=function(war){return(NULL)},
        error=function(err){return(NULL)}
        )
        if(!is.null(fit.out)){
          logliks=c(logliks,fit.out$loglik)
        }else{
          logliks=c(logliks,NA)
        }
      }
      # if some distributions got the same log likelihood, choose one at random
      logliks=logliks+runif(length(logliks),min=0,max=0.001)
      max.index=which(logliks==max(logliks,na.rm=TRUE))
      #print(paste("Selected distribution for taxon",i,"is",distribtypes[max.index]))
      #print(paste("Goodness of fit",logliks[max.index]))
      if(length(max.index) > 0){
        distribs=c(distribs,distribtypes[max.index])
      }else{
        distribs=c(distribs,"NA")
      }
    }else{
      distribs=c(distribs,"NA")
    }
  } # end loop rows
  distribsummary=list(which(distribs=="normal"),which(distribs=="lognormal"),which(distribs=="negative binomial"),which(distribs=="Poisson"),which(distribs=="NA"))
  names(distribsummary)=c("norm","lognorm","negbin","pois","none")
  print(paste("Number of taxa with normal distribution: ",length(distribsummary$norm)))
  print(paste("Number of taxa with lognormal distribution: ",length(distribsummary$lognorm)))
  print(paste("Number of taxa with negative binomial distribution: ",length(distribsummary$negbin)))
  print(paste("Number of taxa with Poisson distribution: ",length(distribsummary$pois)))
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

