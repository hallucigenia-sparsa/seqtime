#' @title Quality scores and plot for estimated interaction matrices
#'
#' @description The predicted time series can be obtained by simulating the time series
#' with Ricker or SOI at once or step-wise by computing the current abundances from the
#' original abundances of the preceding time point. The model parameters are supposed to have
#' been estimated from the time series previously. In case of Ricker, the function estimates carrying capacities
#' as mean abundances. If the original interaction matrix is provided,
#' the mean correlation between estimated and original interaction matrix is computed.
#' Thus, two main quality scores are considered: the mean correlation between predicted and
#' observed time series and, if provided, the mean correlation between the original and the
#' estimated interaction matrix. To avoid explosions in case of Ricker, the estimated interaction matrix can be modified
#' such that all its eigen values are negative using the Schur decomposition.
#'
#' @param oriTS the original time series, with taxa as rows and time points as columns
#' @param A the estimated interaction matrix
#' @param A.ori the original interaction matrix (optional)
#' @param type the model used to predict time series from the interaction matrix, ricker or soi (for soi, predict.stepwise is set to FALSE and noSchur to TRUE)
#' @param m.vector the immigration rates for soi (optional, if not provided set to the proportions in the first sample)
#' @param e.vector the extinction rates for soi (optional, if not provided sampled from the uniform distribution)
#' @param spec.subset either a vector of species indices to keep or a number to indicate how many top-abundant (sum across samples) species should be kept, applied to the matrices and time series
#' @param norm normalize the original time series by dividing each sample by its sum (carried out before filtering species)
#' @param plot plot the number of species versus the mean correlation of predicted and observed time series and the mean auto-correlation
#' @param predict.stepwise if TRUE, the predicted time series is computed step by step, else computed with a call to Ricker
#' @param sim similarity measure to compare predicted and observed time series, either Pearson (r, default) or Kullback-Leibler dissimilarity (kld)
#' @param noSchur do not remove positive eigenvalues
#' @param ignoreExplosion ignore the occurrence of an explosion (for step-wise prediction only)
#' @param sigma noise factor in Ricker
#' @param explosion.bound the explosion boundary in Ricker
#' @return a list with the taxon numbers considered (taxonnum), the mean correlation of observed and predicted time series (meancrosscor), the slope of the former (slope), the value range in the estimated matrix (AestR), the mean correlation between predicted and original matrix (corAAest) and the mean autocorrelations up to lag 5 (meanautocor1 to meanautocor5)
#' @examples
#' N=20
#' A=generateA(N,c=0.1)
#' ts=ricker(N=N,A=A)
#' ts=seqtime::normalize(ts)
#' Aest=limits(ts,verbose=TRUE)
#' out=limitsQuality(ts,A=Aest,A.ori=A, plot=TRUE)
#' @export

limitsQuality<-function(oriTS, A, A.ori=matrix(), type="ricker", m.vector=c(), e.vector=c(), spec.subset=NA, norm=FALSE, plot=FALSE, predict.stepwise=TRUE, sim="r", noSchur=FALSE, ignoreExplosion=FALSE, sigma=-1, explosion.bound=10^8){
  if(ncol(oriTS) < 2){
    stop("Please provide the original time series.")
  }

  if(is.null(A) || ncol(A) < 2){
    stop("Please provide the estimated interaction matrix.")
  }

  if(sim!="kld" && sim!="r"){
    stop("sim should be either r (Pearson) or kld (Kullback-Leibler dissimilarity).")
  }

  # normalize the time series
  if(norm == TRUE){
    oriTS=seqtime::normalize(oriTS)
  }

  # filter out species if requested
  if(!is.na(spec.subset)){
    if(length(spec.subset)==1){
      # only keep top-abundant species
      sorted=sort(apply(oriTS,1,sum),decreasing=TRUE,index.return=TRUE)
      spec.subset=sorted$ix[1:spec.subset]
    }
    if(length(spec.subset)>1){
      oriTS=oriTS[spec.subset,]
      A=A[spec.subset,spec.subset]
      if(ncol(A.ori)>1){
        A.ori=A.ori[spec.subset,spec.subset]
      }
    }
  }

  N=ncol(A)
  corAwithAest=NA
  rangeAest=range(A)
  rangeAest=rangeAest[2]-rangeAest[1]
  slope=NA

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
  # make interaction matrix more robust to explosions (not needed for soc)
  if(noSchur == FALSE && type=="ricker"){
    print("Applying Schur decomposition")
    Amodif=modifyA(A,mode="schur")
  }else{
    Amodif=A
  }
  # loop over the number of species to include (removing the least abundant)
  indicesOK=c()
  percOKlist<-c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0)
  for(percOK in percOKlist){
    indicesOK<-which(rowMeans>percOK*maxAbundance)
    #print(indicesOK)
    if(length(indicesOK)>1){
      specNumberVec<-cbind(specNumberVec,length(indicesOK))
      subTS=oriTS[indicesOK,]
      Amodifsub=Amodif[indicesOK,indicesOK]
      # compute cross-correlation of original and step-wise predicted time series
      if(predict.stepwise == TRUE && type=="ricker"){
        crossres=getCrossCorOriPredStepwise(subTS,A=Amodifsub,sim=sim,sigma=sigma,explosion.bound = explosion.bound, ignoreExplosion=ignoreExplosion)
      }else{
        crossres=getCrossOriPredFull(subTS,A=Amodifsub,sim=sim,sigma=sigma,explosion.bound = explosion.bound, type=type, m.vector=m.vector, e.vector=e.vector)
      }
      crosscorVec=c(crosscorVec,crossres)
      # compute auto-correlations for lag 1 to 5
      if(sim=="kld"){
        autocorrStep1Vec=c(autocorrStep1Vec, getMeanAutoKLD(subTS,lag=1))
        autocorrStep2Vec=c(autocorrStep2Vec, getMeanAutoKLD(subTS,lag=2))
        autocorrStep3Vec=c(autocorrStep3Vec, getMeanAutoKLD(subTS,lag=3))
        autocorrStep4Vec=c(autocorrStep4Vec, getMeanAutoKLD(subTS,lag=4))
        autocorrStep5Vec=c(autocorrStep5Vec, getMeanAutoKLD(subTS,lag=5))
      }else if(sim=="r"){
        autocorrStep1Vec=c(autocorrStep1Vec, getMeanAutocor(subTS,lag=1))
        autocorrStep2Vec=c(autocorrStep2Vec, getMeanAutocor(subTS,lag=2))
        autocorrStep3Vec=c(autocorrStep3Vec, getMeanAutocor(subTS,lag=3))
        autocorrStep4Vec=c(autocorrStep4Vec, getMeanAutocor(subTS,lag=4))
        autocorrStep5Vec=c(autocorrStep5Vec, getMeanAutocor(subTS,lag=5))
      }
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

  # compute slope between species number and mean cross-correlation
  reg.data=data.frame(specNumberVec,crosscorVec)
  linreg = lm(formula = crosscorVec~specNumberVec)
  slope=linreg$coefficients[2]

  # compute slope between species number and mean auto-correlation of lag 1
  if(length(which(is.na(autocorrStep1Vec)))<length(autocorrStep1Vec)){
    auto.reg.data=data.frame(specNumberVec,autocorrStep1Vec)
    auto.linreg = lm(formula = autocorrStep1Vec~specNumberVec)
    autoslope=auto.linreg$coefficients[2]
  }else{
    autoslope=NA
  }

  # prepare result object
  res=list(specNumberVec,crosscorVec,autocorrStep1Vec,autocorrStep2Vec,autocorrStep3Vec,autocorrStep4Vec,autocorrStep5Vec)
  names(res)=c("taxonnum","meancrosscor","meanautocor1","meanautocor2","meanautocor3","meanautocor4","meanautocor5")

  if(plot == TRUE){
    # do the quality plot
    mat=matrix(unlist(res),nrow=length(specNumberVec),ncol=length(res))
    # exclude taxon number column
    mat=mat[,2:ncol(mat)]

    main="Quality of estimated interaction matrix"
    colnames(mat)=c("cor(ori.ts,pred.ts)","autocor(ori.ts) lag 1","autocor(ori.ts) lag 2","autocor(ori.ts) lag 3","autocor(ori.ts) lag 4","autocor(ori.ts) lag 5")
    my.colors = c("red","blue","blue4","violet","cyan","cadetblue")

    # the extra 50 is place for the legend
    plot(specNumberVec,mat[,ncol(mat)],xlim=c(0,N+100),ylim = c(-1,1), xlab = "Number of taxa (in decreasing abundance)", ylab = "Mean correlation", main = main, type = "o", col = my.colors[ncol(mat)])
    # plot the cross-correlation prediction last
    for(i in (ncol(mat)-1):1){
      lines(specNumberVec,mat[,i], col = my.colors[i], type="o")
    }
    legend(x="right",colnames(mat), lty = rep(1,ncol(mat)), col = my.colors, merge = TRUE, bg = "white", text.col="black")
  }
  res["corAAest"]=corAwithAest
  res["AestR"]=rangeAest
  res["slope"]=slope
  res["autoslope"]=autoslope
  return(res)
}

# Generate predicted time series by calling Ricker or SOC with the given interaction matrix.
#
# oriTS: original time series, taxa are rows and time points are columns
# A: estimated interaction matrix
# sim: similarity to use to assess discrepancy between observed and predicted time series
# sigma: the noise term
# explosion.bound: the explosion boundary of Ricker
# type: ricker or soc model
# m.vector: the vector of immigration probabilities for SOC (optional)
# e.vector: the vector of extinction probabilities for SOC (optional)
#
# returns the mean cross-correlation of the predicted time series
getCrossOriPredFull<-function(oriTS, A, sim="r", sigma=-1, explosion.bound=10^8, type="ricker", m.vector=c(), e.vector=c()){
  y=oriTS[,1]
  tend=ncol(oriTS)
  N=nrow(oriTS)
  if(type=="ricker"){
    # estimate carrying capacity as the mean of the time series
    K=apply(oriTS,1,mean)
    stable=testStability(A, method="ricker", K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
    if(stable){
      predTS=ricker(N=N, A=A, sigma=sigma, y=y, tend=tend, K=K, explosion.bound=explosion.bound)
    }else{
      stop("Cannot generate predicted time series because of explosion in Ricker.")
    }
  }else{
    sums=apply(oriTS,2,sum)
    I=round(mean(sums))
    if(length(m.vector)==0){
      m.vector=oriTS[,1]/sum(oriTS[,1])
    }
    if(length(e.vector)==0){
      e.vector=runif(N,min=0,max=1)
    }
    predTS=soc(N, I, A, m.vector=m.vector, e.vector=e.vector, tend)
  }
  if(sim=="kld"){
    crossCor<-getMeanKLD(oriTS[,2:tend],predTS[,2:tend],lag=0)
  }else if(sim=="r"){
    crossCor<-getMeanCrosscor(oriTS[,2:tend],predTS[,2:tend],lag=0)
  }
  return(crossCor)
}

# Generate predicted time series step-wise by calculating the abundances for
# each time point from the original abundances of the preceding time point
# using Ricker with the given interaction matrix.
#
# oriTS: original time series, taxa are rows and time points are columns
# A: estimated interaction matrix
# sim: similarity to use to assess discrepancy between observed and predicted time series
# lag: the predecing time point to consider (1 = one step before current, 2 = two steps before current etc.)
# sigma: the noise term
# explosion.bound the explosion boundary of Ricker
#
# returns the mean cross-correlation of the predicted time series
getCrossCorOriPredStepwise<-function(oriTS, A, sim="r", lag=1,sigma=-1, explosion.bound=10^8, ignoreExplosion=FALSE){
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
    if(max(y) > explosion.bound){
      if(ignoreExplosion == FALSE){
        stop(paste("Explosion for taxon",which(y==max(y)),"in step-wise prediction."))
      }else{
        # even if ignored, report explosion
        warning(paste("Explosion for taxon",which(y==max(y)),"in step-wise prediction."))
      }
    }
    predTS[,t]=y
  }
  if(sim=="kld"){
    crossCor<-getMeanKLD(oriTS[,2:tend],predTS[,2:tend],lag=0)
  }else if(sim=="r"){
    crossCor<-getMeanCrosscor(oriTS[,2:tend],predTS[,2:tend],lag=0)
  }
  return(crossCor)
}

# Get the mean KL dissimilarity between the
# time series, with a shift given by lag.
# Taxa are rows and time points columns.
getMeanAutoKLD<-function(ts, lag=1){
  tend=ncol(ts)
  # correlation between matrices is computed column-wise by default, so transpose
  return(getMeanKLD(ts[,1:(tend-lag)],ts[,(1+lag):tend]))
}

# Get the mean KL dissimilarity between the observed and predicted
# time series, with a shift given by lag.
# oriTS, predTS: observed and predicted community time series with taxa as rows and time points as columns
# lag: lag between observed and predicted
getMeanKLD<-function(oriTS, predTS, lag=0){
  tend=ncol(oriTS)
  # compute probabilities
  probabOri=getProbab(oriTS[,1:(tend-lag)])
  probabPred=getProbab(predTS[,(1+lag):tend])
  klds=c()
  for(i in 1:nrow(oriTS)){
    kld = get.kld(oriTS[i,],predTS[i,],probab=TRUE)
    klds=c(klds,kld)
  }
  return(mean(klds))
}

##############################################################################
# Compute Kullback Leibler dissimilarity (symmetrized divergence) between two vectors.
#
# Argument:
# x = vector with non-negative numbers
# y = vector with non-negative numbers
# pseudocount = this value is added to each zero
# probab = if false, x and y are divided by their sum so that they sum to one
#
# Value:
# Kullbakc Leibler dissimilarity
#
# Note:
# Equation: D(x,y) = SUM(x_i*log(x_i/y_i) + y_i*log(y_i/x_i))
# taken from "Caution! Compositions! Can constraints on omics data lead analyses astray?"
# David Lovell et al., Report Number EP10994
#
##############################################################################
get.kld=function(x,y, pseudocount=0.00000001, probab=FALSE){
  if(length(x) != length(y)){
    stop("The two vectors should have the same length!")
  }
  x[x==0]=pseudocount
  y[y==0]=pseudocount
  dis = 0
  if(probab==FALSE){
    x = x/sum(x)
    y = y/sum(y)
  }
  for(i in 1:length(x)){
    if(!is.nan(x[i]) && !is.nan(y[i])){
      ratioxy = log(x[i]/y[i])
      ratioyx = log(y[i]/x[i])
      dis = x[i]*ratioxy+y[i]*ratioyx + dis
    }
  }
  dis
}

# TS is a matrix with taxa as rows and time points as columns
# Each row is divided by its sum.
getProbab<-function(TS){
  rowSums=apply(TS,1,sum)
  return(TS/rowSums)
}

# Compute mean Cameron & Windmeijer R2 across all time series.
# Cameron & Windmeijer (1997) An R-squared measure of goodness of fit for some common nonlinear regression models.
# Journal of Econometrics 77, pp. 329-342.
#
# oriTS, predTS: observed and predicted community time series with taxa as rows and time points as columns
# lag: lag between observed and predicted
getMeanCWR2<-function(oriTS, predTS, lag=0){
  tend=ncol(oriTS)
  # compute probabilities
  probabOri=getProbab(oriTS[,1:(tend-lag)])
  probabPred=getProbab(predTS[,(1+lag):tend])
  cwr2s=c()
  mean=mean(apply(probabOri,1,mean))
  meanVec=rep(mean,ncol(oriTS))
  for(i in 1:nrow(oriTS)){
    cwr2= 1 - (gKL(probabOri[i,], probabPred[i,])/gKL(probabOri[i,], meanVec))
    cwr2s=c(cwr2s,cwr2)
  }
  return(mean(cwr2s))
}

# Compute generalized Kullback-Leibler divergence for two vectors representing probability distributions.
# Cameron & Windmeijer (1997) An R-squared measure of goodness of fit for some common nonlinear regression models.
# Journal of Econometrics 77, pp. 329-342.
gKL<-function(vec1,vec2){
  ratios=vec1/vec2
  # log: by default natural
  return(2*sum(log(ratios)))
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
  # correlation between matrices is computed column-wise by default, so transpose
  Rcross=cor(t(oriTS[,1:(tend-lag)]),t(predTS[,(1+lag):tend]))
  return(mean(diag(Rcross), na.rm=TRUE))
}

# Get the mean KL divergence between the
# time series, with a shift given by lag.
# Taxa are rows and time points columns.
getMeanAutoCWR2<-function(ts, lag=1){
  tend=ncol(ts)
  # correlation between matrices is computed column-wise by default, so transpose
  return(getMeanCWR2(ts[,1:(tend-lag)],ts[,(1+lag):tend]))
}

# Get the mean auto-correlation between the
# time series, with a shift given by lag.
# Taxa are rows and time points columns.
getMeanAutocor<-function(ts,lag=1){
  tend=ncol(ts)
  # correlation between matrices is computed column-wise by default, so transpose
  Rauto=cor(t(ts[,1:(tend-lag)]),t(ts[,(1+lag):tend]))
  return(mean(diag(Rauto)))
}
