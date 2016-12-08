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
