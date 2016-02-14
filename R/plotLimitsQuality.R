#' Quality plot for interaction matrices estimated with LIMITS
#'
#' To avoid explosions, the estimated interaction matrix is modified such
#' that all its eigen values are negative using the Schur decomposition. The predicted time series is
#' obtained step-wise, by computing the current abundances from the
#' original abundances of the preceding time point and the given
#' estimated interaction matrix with (by default noise-free) Ricker. Carrying capacities are
#' estimated as mean abundances. If the original interaction matrix is provided,
#' the mean correlation between estimated and original interaction matrix is computed.
#'
#' @param oriTS the original time series, with taxa as rows and time points as columns
#' @param A the estimated interaction matrix
#' @param A.ori the original interaction matrix (optional)
#' @param predict.stepwise if TRUE, the predicted time series is computed step by step, else computed with a call to Ricker
#' @param noSchur do not remove positive eigenvalues
#' @param ignoreExplosion ignore the occurrence of an explosion (for step-wise prediction only)
#' @param sigma noise factor in Ricker
#' @param explosion.bound the explosion boundary in Ricker
#' @return a list with the taxon numbers considered (taxonnum), the mean correlation of observed and predicted time series (meancrosscor) and the mean autocorrelations up to lag 5 (meanautocor1 to meanautocor5)
#' @examples
#' N=20
#' A=generateA(N,c=0.1)
#' ts=ricker(N=N,A=A)
#' Aest=limits(ts,verbose=TRUE)
#' out=plotLimitsQuality(ts,A=Aest,A.ori=A)
#' @export

plotLimitsQuality<-function(oriTS, A, A.ori=matrix(), predict.stepwise=TRUE, noSchur=FALSE, ignoreExplosion=FALSE, sigma=-1, explosion.bound=10^8){
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
  autocorrStep1Vec=c()
  autocorrStep2Vec=c()
  autocorrStep3Vec=c()
  autocorrStep4Vec=c()
  autocorrStep5Vec=c()
  # remove problematic species from the interaction matrix
  if(noSchur == FALSE){
    print("Applying Schur decomposition")
    Amodif=modifyA(A,mode="schur")
  }else{
    Amodif=A
  }
  # loop over the number of species to include (removing the least abundant)
  percOKlist<-c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0)
  for(percOK in percOKlist){
    indicesOK<-which(rowMeans>percOK*maxAbundance)
    if(length(indicesOK)>1){
      specNumberVec<-cbind(specNumberVec,length(indicesOK))
      subTS=oriTS[indicesOK,]
      Amodifsub=Amodif[indicesOK,indicesOK]
      # compute cross-correlation of original and step-wise predicted time series
      if(predict.stepwise == TRUE){
        crossres=getCrossCorOriPredStepwise(subTS,A=Amodifsub,sigma=sigma,explosion.bound = explosion.bound, ignoreExplosion=ignoreExplosion)
      }else{
        crossres=getCrossOriPredFull(subTS,A=Amodifsub,sigma=sigma,explosion.bound = explosion.bound)
      }
      crosscorVec=c(crosscorVec,crossres)
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
  autocorrStep1Vec=autocorrStep1Vec[-1]
  autocorrStep2Vec=autocorrStep2Vec[-1]
  autocorrStep3Vec=autocorrStep3Vec[-1]
  autocorrStep4Vec=autocorrStep4Vec[-1]
  autocorrStep5Vec=autocorrStep5Vec[-1]

  # prepare result object
  res=list(specNumberVec,crosscorVec,autocorrStep1Vec,autocorrStep2Vec,autocorrStep3Vec,autocorrStep4Vec,autocorrStep5Vec)
  names(res)=c("taxonnum","meancrosscor","meanautocor1","meanautocor2","meanautocor3","meanautocor4","meanautocor5")

  # do the quality plot
  mat=matrix(unlist(res),nrow=length(specNumberVec),ncol=length(res))
  # exclude taxon number column
  mat=mat[,2:ncol(mat)]
  colnames(mat)=c("cor(ori.ts,pred.ts)","autocor(ori.ts) lag 1","autocor(ori.ts) lag 2","autocor(ori.ts) lag 3","autocor(ori.ts) lag 4","autocor(ori.ts) lag 5")

  my.colors = c("red","blue","blue4","violet","cyan","cadetblue")
  # the extra 50 is place for the legend
  plot(specNumberVec,mat[,1],xlim=c(0,N+100),ylim = c(-1,1), xlab = "Number of taxa considered", ylab = "Mean correlation", main = "Quality of estimated interaction matrix", type = "o", col = my.colors[1])
  for(i in 2:ncol(mat)){
    lines(specNumberVec,mat[,i], col = my.colors[i], type="o")
  }
  legend(x="right",colnames(mat), lty = rep(1,ncol(mat)), col = my.colors, merge = TRUE, bg = "white", text.col="black")

  return(res)
}

# Generate predicted time series by calling Ricker with the given interaction matrix.
#
# oriTS: original time series, taxa are rows and time points are columns
# A: estimated interaction matrix
# sigma: the noise term
# explosion.bound the explosion boundary of Ricker
#
# returns the mean cross-correlation of the predicted time series
getCrossOriPredFull<-function(oriTS, A, sigma=-1, explosion.bound=10^8){
  y=oriTS[,1]
  tend=ncol(oriTS)
  N=nrow(oriTS)
  # estimate carrying capacity as the mean of the time series
  K=apply(oriTS,1,mean)
  stable=testStability(A, method="ricker", K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
  if(stable){
    predTS=ricker(N=N, A=A, sigma=sigma, y=y, tend=tend, K=K, explosion.bound=explosion.bound)
    crossCor<-getMeanCrosscor(oriTS[,2:tend],predTS[,2:tend],lag=0)
    return(crossCor)
  }else{
    stop("Cannot generate predicted time series because of explosion in Ricker.")
  }
}

# Generate predicted time series step-wise by calculating the abundances for
# each time point from the original abundances of the preceding time point
# using Ricker with the given interaction matrix.
#
# oriTS: original time series, taxa are rows and time points are columns
# A: estimated interaction matrix
# lag: the predecing time point to consider (1 = one step before current, 2 = two steps before current etc.)
# sigma: the noise term
# explosion.bound the explosion boundary of Ricker
#
# returns the mean cross-correlation of the predicted time series
getCrossCorOriPredStepwise<-function(oriTS,A,lag=1,sigma=-1, explosion.bound=10^8, ignoreExplosion=FALSE){
  tend=ncol(oriTS)
  N=nrow(oriTS)
  # estimate carrying capacity as the mean of the time series
  K=apply(oriTS,1,mean)
  predTS = matrix(0,nrow=N,ncol=tend)
  predTS[,1]=as.numeric(oriTS[,1])
  # generate predicted time series from original time series and estimated interaction matrix step by step
  for(t in 2:tend){
    prevT<- t-lag
    # noise factor is optionally considered
    b=rep(1,N)
    if(sigma > 0){
      b=rlnorm(N,meanlog=0,sdlog=sigma)
    }
    y<-b*as.numeric(oriTS[,prevT])*exp(A%*%(as.numeric(oriTS[,prevT])-K))
    if(max(y) > explosion.bound && ignoreExplosion == FALSE){
      # report which species explodes
      stop(paste("Explosion for taxon",which(y==max(y)),"in step-wise prediction."))
    }
    predTS[,t]=y
  }
  crossCor<-getMeanCrosscor(oriTS[,2:tend],predTS[,2:tend],lag=0)
  return(crossCor)
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
