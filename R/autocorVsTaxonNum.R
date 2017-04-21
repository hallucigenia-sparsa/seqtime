#' @title Compute slope for increasing number of taxa versus their autocorrelation
#' @description The taxa are added in decreasing order of abundance.
#'
#' @param x a matrix with objects as rows and time points as columns
#' @param lag the lag for which autocorrelation is computed
#' @param plot plot the increasing number of taxa versus their autocorrelation
#' @return a list with the taxon number, autocorrelations, slope and p-value
#'
#' @examples
#' N=10
#' A=generateA(N,c=0.1)
#' soi.ts=soi(N=N,I=1000,A=A,tend=100)
#' res=autocorVsTaxonNum(soi.ts, plot=TRUE)
#' @export

autocorVsTaxonNum<-function(ts, lag=1, plot=FALSE){
  tend=ncol(ts)
  rowMeans=apply(ts,1,mean)
  sortedTaxa=sort(rowMeans,decreasing=TRUE,index.return=TRUE)$ix
  indicesOK=c()
  autocorrels=c()
  specNumberVec=c()
  for(taxon in sortedTaxa){
    indicesOK=c(indicesOK, taxon)
    if(length(indicesOK)>1){
      specNumberVec<-c(specNumberVec,length(indicesOK))
      #print(indicesOK)
      subTS=ts[indicesOK,]
      # correlation between matrices is computed column-wise by default, so transpose
      Rauto=cor(t(subTS[,1:(tend-lag)]),t(subTS[,(1+lag):tend]))
      autocorrels=c(autocorrels, mean(diag(Rauto)))
    } # got more than one species
  } # loop over species numbers to consider

  # compute slope between species number and mean auto-correlation
  if(length(which(is.na(autocorrels)))<length(autocorrels)){
    auto.reg.data=data.frame(specNumberVec,autocorrels)
    auto.linreg = lm(formula = autocorrels~specNumberVec, data=auto.reg.data)
    sum=summary(auto.linreg)
    pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
    autoslope=auto.linreg$coefficients[2]
  }else{
    pval=NA
    autoslope=NA
  }

  if(plot==TRUE){
    plot(specNumberVec,autocorrels,main=paste("Mean autocorrelation for lag ",lag," versus taxon number\nSlope: ",round(autoslope,4),", p-value: ",round(pval,4),sep=""), xlab="Number of taxa added in decreasing order of abundance", ylab="Mean autocorrelation", ylim=c(0,1))
    if(!is.na(autoslope)){
      abline(auto.linreg,bty="n",col="red")
    }
  }

  res=list(specNumberVec,autocorrels,autoslope,pval)
  names(res)=c("taxonNumber","autocor","slope","pval")
  return(res)
}
