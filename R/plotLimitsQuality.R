#' Quality plot for interaction matrices estimated with LIMITS
#'
#' To avoid explosions, the estimated interaction matrix is modified such
#' that all its eigen values are negative. The predicted time series is
#' obtained step-wise, by computing the current abundances from the
#' original abundances of the preceding time point and the given
#' estimated interaction matrix with (noise-free) Ricker. Carrying capacities are
#' estimated as mean abundances. If the original interaction matrix is provided,
#' the mean correlation between estimated and original interaction matrix is computed.
#'
#' @param oriTS the original time series, with taxa as rows and time points as columns
#' @param A the estimated interaction matrix
#' @param A.ori the original interaction matrix (optional)
#' @param sigma noise term for time series generation with ricker
#' @param explosion.bound explosion boundary for time series generation with ricker
#' @return a list with the taxon numbers considered (taxonnum), the mean correlation of observed and predicted time series from lag 0 to lag 4 (meancrosscor and meancrosscor1 to meancrosscor4) and the mean autocorrelations up to lag 4 (meanautocor1 to meanautocor4)
#' @export

plotLimitsQuality<-function(oriTS, A, A.ori=matrix()){
  if(ncol(oriTS) < 2){
    stop("Please provide the original time series.")
  }
  if(ncol(A) < 2){
    stop("Please provide the estimated interaction matrix.")
  }
  N=ncol(A)
  if(ncol(A.ori) > 1){
    # compute correlations between the same columns only
    corAwithAest=sum(diag(cor(A,A.ori)))/N
    print(paste("Average correlation between columns of the original and the estimated interaction matrix:",corAwithAest))
  }
  # prepare vectors
  rowMeans=apply(oriTS,1,mean)
  maxAbundance=max(rowMeans)
  tend=ncol(oriTS)
  specNumberVec=c()
  crosscorVec=c()
  crosscorStep1Vec=c()
  crosscorStep2Vec=c()
  crosscorStep3Vec=c()
  crosscorStep4Vec=c()
  autocorrStep1Vec=c()
  autocorrStep2Vec=c()
  autocorrStep3Vec=c()
  autocorrStep4Vec=c()
  autocorrStep5Vec=c()
  # remove problematic species from the interaction matrix
  Amodif=modifyA(A,mode="schur")
  # loop over the number of species to include in the order of their abundance
  percOKlist<-c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0)
  for(percOK in percOKlist){
    indicesOK<-which(rowMeans>percOK*maxAbundance)
    if(length(indicesOK)>1){
      specNumberVec<-cbind(specNumberVec,length(indicesOK))
      subTS=oriTS[indicesOK,]
      Amodifsub=Amodif[indicesOK,indicesOK]
      # compute cross-correlation of original and step-wise predicted time series
      res=getCrossCorOriPred(subTS,A=Amodifsub)
      crosscorVec=c(crosscorVec,res$meancrosscor)
      crosscorStep1Vec=c(crosscorStep1Vec,res$meancrosscor1)
      crosscorStep2Vec=c(crosscorStep2Vec,res$meancrosscor2)
      crosscorStep3Vec=c(crosscorStep3Vec,res$meancrosscor3)
      crosscorStep4Vec=c(crosscorStep4Vec,res$meancrosscor4)
      # compute auto-correlations for lag 1 to 5
      autocorrStep1Vec=c(autocorrStep1Vec, getMeanAutocor(subTS,lag=1))
      autocorrStep2Vec=c(autocorrStep2Vec, getMeanAutocor(subTS,lag=2))
      autocorrStep3Vec=c(autocorrStep3Vec, getMeanAutocor(subTS,lag=3))
      autocorrStep4Vec=c(autocorrStep4Vec, getMeanAutocor(subTS,lag=4))
      autocorrStep5Vec=c(autocorrStep5Vec, getMeanAutocor(subTS,lag=5))
    } # got more than one species
  } # loop over species numbers to consider

  # remove first element
  specNumberVec=specNumberVec[-1]
  crosscorVec=crosscorVec[-1]
  crosscorStep1Vec=crosscorStep1Vec[-1]
  crosscorStep2Vec=crosscorStep2Vec[-1]
  crosscorStep3Vec=crosscorStep3Vec[-1]
  crosscorStep4Vec=crosscorStep4Vec[-1]
  autocorrStep1Vec=autocorrStep1Vec[-1]
  autocorrStep2Vec=autocorrStep2Vec[-1]
  autocorrStep3Vec=autocorrStep3Vec[-1]
  autocorrStep4Vec=autocorrStep4Vec[-1]
  autocorrStep5Vec=autocorrStep5Vec[-1]

  # prepare result object
  res=list(specNumberVec,crosscorVec,crosscorStep1Vec,crosscorStep2Vec,crosscorStep3Vec,crosscorStep4Vec,autocorrStep1Vec,autocorrStep2Vec,autocorrStep3Vec,autocorrStep4Vec)
  names(res)=c("taxonnum","meancrosscor","meancrosscor1","meancrosscor2","meancrosscor3","meancrosscor4","meanautocor1","meanautocor2","meanautocor3","meanautocor4")

  # do the quality plot
  mat=matrix(unlist(res),nrow=length(specNumberVec),ncol=length(res))
  # exclude taxon number column
  mat=mat[,2:ncol(mat)]
  colnames(mat)=c("cor(ori.ts,pred.ts)", "cor(ori.ts,pred.ts) 1 step","cor(ori.ts,pred.ts) 2 steps", "cor(ori.ts,pred.ts) 3 steps","cor(ori.ts,pred.ts) 4 steps","autocor(ori.ts) lag 1","autocor(ori.ts) lag 2","autocor(ori.ts) lag 3","autocor(ori.ts) lag 4")

  my.colors = c("red","orange","burlywood","brown","chocolate","blue","blue4","violet","cyan")
  # the extra 50 is place for the legend
  plot(specNumberVec,mat[,1],xlim=c(0,N+100),ylim = c(-1,1), xlab = "Number of taxa considered", ylab = "Mean correlation", main = "Quality of estimated interaction matrix", type = "o", col = my.colors[1])
  for(i in 2:ncol(mat)){
    lines(specNumberVec,mat[,i], col = my.colors[i], type="o")
  }
  legend(x="right",colnames(mat), lty = rep(1,ncol(mat)), col = my.colors, merge = TRUE, bg = "white", text.col="black")

  return(res)
}

# Generate predicted time series step-wise by calculating the abundances for
# each time point from the original abundances of the preceding time point
# using Ricker with the given interaction matrix.
# In addition, additional time series are created by predicting
# the next step from the predicted time series, for up to 4 steps.
#
# oriTS: original time series, taxa are rows and time points are columns
# A: estimated interaction matrix
#
# a list with the mean cross-correlation of the predicted time series
getCrossCorOriPred<-function(oriTS,A,lag=1){
  tend=ncol(oriTS)
  # estimate carrying capacity as the mean of the time series
  K=apply(oriTS,1,mean)
  oriTS=t(oriTS)
  predTS <- oriTS[1,]
  predTS2 <- t(oriTS[1:2,])
  predTS3 <- t(oriTS[1:3,])
  predTS4 <- t(oriTS[1:4,])
  predTS5 <- t(oriTS[1:5,])
  # generate predicted time series from original time series and estimated interaction matrix step by step
  for(t in 2:tend){
    prevT<- t-lag
    predTS<-cbind(predTS,oriTS[prevT,]*exp(A%*%(oriTS[prevT,]-K)))
    predTS2<-cbind(predTS2,predTS[,prevT]*exp(A%*%(predTS[,prevT]-K)))
    predTS3<-cbind(predTS3,predTS2[,prevT]*exp(A%*%(predTS2[,prevT]-K)))
    predTS4<-cbind(predTS4,predTS3[,prevT]*exp(A%*%(predTS3[,prevT]-K)))
    predTS5<-cbind(predTS5,predTS4[,prevT]*exp(A%*%(predTS4[,prevT]-K)))
  }
  predTS<-t(predTS)
  predTS2<-t(predTS2[,1:tend])
  predTS3<-t(predTS3[,1:tend])
  predTS4<-t(predTS4[,1:tend])
  predTS5<-t(predTS5[,1:tend])
  crossCor1<-getMeanCrosscor(oriTS[2:tend,],predTS[2:tend,],lag=0)
  crossCor2<-getMeanCrosscor(oriTS[3:tend,],predTS[3:tend,],lag=0)
  crossCor3<-getMeanCrosscor(oriTS[4:tend,],predTS[4:tend,],lag=0)
  crossCor4<-getMeanCrosscor(oriTS[5:tend,],predTS[5:tend,],lag=0)
  crossCor5<-getMeanCrosscor(oriTS[6:tend,],predTS[6:tend,],lag=0)
  res=list(crossCor1,crossCor2,crossCor3,crossCor4,crossCor5)
  names(res)=c("meancrosscor","meancrosscor1","meancrosscor2","meancrosscor3","meancrosscor4")
  return(res)
}

# Get the mean auto-correlation between the
# time series, with a shift given by lag.
# Taxa are rows and time points columns.
# If no shift is needed, set lag to 0.
getMeanCrosscor<-function(oriTS,predTS,lag=1){
  if(ncol(oriTS) != ncol(predTS)){
    stop("Original and predicted time series do not have the same number of time points!")
  }
  tend=ncol(oriTS)
  # correlation between matrices is computed column-wise
  Rcross=cor(t(oriTS[,1:(tend-lag)]),t(predTS[,(1+lag):tend]))
  return(mean(diag(Rcross)))
}

# Get the mean auto-correlation between the
# time series, with a shift given by lag.
# Taxa are rows and time points columns.
getMeanAutocor<-function(ts,lag=1){
  tend=ncol(ts)
  # correlation between matrices is computed column-wise
  Rauto=cor(t(ts[,1:(tend-lag)]),t(ts[,(1+lag):tend]))
  return(mean(diag(Rauto)))
}
