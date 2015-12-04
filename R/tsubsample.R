#' Subsample a time series
#'
#' Subsample the given time series down to a target time point number using equal intervals.
#' Assumes that time points are columns and time points are equidistant.
#' @param x a matrix
#' @param num the resulting number of time points (defaults to 100)
#' @export

tsubsample<-function(x, num=100){
  timepoints=ncol(x)
  intervals=round(timepoints/num)
  samples=seq(1,timepoints,by=intervals)
  return(x[,samples])
}
