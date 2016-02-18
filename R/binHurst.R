#' Given a community time series, bin taxa according to their Hurst exponent
#'
#' Community taxa can be placed into a low, middle, high and very high Hurst category. The
#' thresholds between categories can be user-defined.
#'
#' @param the time series with taxa as rows and samples as columns
#' @return the indices of the taxa in each Hurst category, labeled lowH, middleH, highH and veryhighH
#' @examples
#' out=binHurst(simuntb(N=100))
#' length(out$veryhighH)
#' length(out$highH)
#' length(out$middleH)
#' length(out$lowH)
#' @export

binHurst<-function(x, thresholds=c(0.5,0.7,0.9)){
  lowH=c()
  middleH=c()
  highH=c()
  veryHighH=c()
  # loop over taxa in time series and compute their maximum autocorrelation for lags > 0
  for(i in 1:nrow(x)){
    h=FGN::HurstK(as.numeric(x[i,]))
    if(!is.na(h)){
      if(h < thresholds[1]){
        lowH=c(lowH,i)
      }
      else if(h >= thresholds[1] && h < thresholds[2]){
        middleH=c(middleH,i)
      }else if(h >= thresholds[2] && h < thresholds[3]){
        highH=c(highH,i)
      }else if(h >= thresholds[3]){
        veryHighH=c(veryHighH,i)
      }
    } # skip NA
  } # loop taxa
  res=list(lowH,middleH,highH,veryHighH)
  names(res)=c("lowH","middleH","highH","veryhighH")
  return(res)
}
