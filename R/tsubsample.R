#' @title Subsample Time Series
#' @description Subsample a time series.
#' @details Subsample the given time series down to a target time point number using equal intervals or select the time points with the given interval. Assumes that time points are columns and time points are equidistant.
#' @param x a matrix
#' @param target target
#' @param interval interval
#' @export
tsubsample <- function(x, target=NA, interval=NA){
  timepoints=ncol(x)
  if(!is.na(target)){
    interval=round(timepoints/target)
  }else if(is.na(interval)){
    stop("Either provide a target or an interval.")
  }
  samples=seq(1,timepoints,by=interval)
  return(x[,samples])
}
