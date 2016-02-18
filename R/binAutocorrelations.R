#' Bin taxa by maximal autocorrelation
#'
#' The taxa are binned according to their maximal autocorrelation into
#' 4 categories with low, middle, high and very high autocorrelation.
#' The autocorrelation at lag = 0 is ignored. The binning thresholds
#' can be altered.
#'
#' @param x the time series with taxa as rows and samples as columns
#' @param thresholds binning thresholds
#' @return a list with the indices of taxa in each autocorrelation bin
#' @examples
#' out=binAutocorrelations(simuntb(N=100))
#' length(out$lowmaxautocor)
#' length(out$highmaxautocor)
#' length(out$veryhighmaxautocor)
#' @export

binAutocorrelations<-function(x, thresholds=c(0.3, 0.5, 0.8)){
  veryhighmaxautocor=c()
  highmaxautocor=c()
  middlemaxautocor=c()
  lowmaxautocor=c()
  # loop over taxa in time series and compute their maximum autocorrelation for lags > 0
  for(i in 1:nrow(x)){
    acfRes=stats::acf(as.numeric(x[i,]),plot=FALSE)
    # skip perfect autocorrelation at lag 0
    maxauto=max(acfRes$acf[2:length(acfRes$acf)])
    if(!is.na(maxauto)){
      if(maxauto < thresholds[1]){
        lowmaxautocor=c(lowmaxautocor,i)
      }else if(maxauto >= thresholds[1] && maxauto < thresholds[2]){
        middlemaxautocor=c(middlemaxautocor,i)
      }else if(maxauto >= thresholds[2] && maxauto < thresholds[3]){
        highmaxautocor=c(highmaxautocor,i)
      }else if(maxauto >= thresholds[3]){
        veryhighmaxautocor=c(veryhighmaxautocor,i)
      }
    } # skip NA
  } # loop taxa
  res=list(lowmaxautocor,middlemaxautocor,highmaxautocor,veryhighmaxautocor)
  names(res)=c("lowmaxautocor","middlemaxautocor","highmaxautocor","veryhighmaxautocor")
  return(res)
}
