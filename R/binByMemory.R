#' @title Bin rows in a matrix by memory.
#'
#' @description Memory is assessed either with the Hurst exponent or with the maximum auto-
#' correlation (omitting lag zero). Hurst exponents are computed with
#' function hurstexp in package pracma. The simple R/S Hurst exponent is selected.
#'
#' @details The function returns an object that groups matrix row indices in different levels.
#' Rows that could not be binned are classified in bin na.
#'
#' @param x a matrix with objects as rows and time points as columns
#' @param thresholds the thresholds by which to bin
#' @param method autocor or hurst, computing the maximum auto-correlation (omitting lag=0) and the Hurst exponent, respectively
#' @param d window size for Hurst exponent computation with function hurstexp
#' @return list with row indices assigned to bins
#' @examples
#' N=10
#' ricker.out=ricker(N,generateA(N),K=rep(0.01,N))
#' memBinned=binByMemory(ricker.out)
#' memBinned.bars=c(length(memBinned$autocor0.8Inf),length(memBinned$autocor0.50.8))
#' names(memBinned.bars)=names(memBinned)
#' barplot(memBinned.bars,main="Auto-correlation bins", ylab="Number of rows")
#' @export

binByMemory<-function(x, thresholds=c(0.3, 0.5, 0.8), method="autocor", d=round(ncol(x)/3)){
  thresholds=c(-Inf,thresholds,Inf)
  res=list()
  # loop rows in x
  for(i in 1:nrow(x)){
    if(method=="hurst"){
      value=pracma::hurstexp(as.numeric(x[i,]),d=d, display=FALSE)
      value=value$Hs
    }else if(method=="autocor"){
      acfRes=stats::acf(as.numeric(x[i,]),plot=FALSE)
      # skip perfect autocorrelation at lag 0
      value=max(acfRes$acf[2:length(acfRes$acf)])
    }else{
      stop("Method ",method," not supported.")
    }
    if(!is.na(value)){
      # loop thresholds
      for(tindex in (1:(length(thresholds)-1))){
        if(value >= thresholds[tindex] && value < thresholds[(tindex+1)]){
          lowerBound=thresholds[tindex]
          if(tindex==1){
            lowerBound="negInf"
          }
          name=paste(method,lowerBound,thresholds[(tindex+1)],sep="")
          if(name %in% names(res)){
            res[[name]]=c(res[[name]],i)
          }else{
            res[[name]]=c(i)
          }
        }
      }
    } # skip NA
    else{
      name="na"
      if(name %in% names(res)){
        res[[name]]=c(res[[name]],i)
      }else{
        res[[name]]=c(i)
      }
    }
  } # loop taxa
  return(res)
}
