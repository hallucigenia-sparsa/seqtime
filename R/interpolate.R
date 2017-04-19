#' @title Time Series Interpolation
#' @description Wrapper function to interpolate a time series.
#'
#' @details If no interval is provided, the intervals in the time vector are computed and the most frequent one is taken as the interval.
#' Note that natural, fmm, parabola and periodic can introduce negative values and are therefore not recommended for taxon abundance data. The default is stineman.
#' @param x the time series matrix, rows are objects and columns are time points
#' @param interval the target intervals to be present after interpolation
#' @param time.index the row index holding time points
#' @param time.vector the vector holding time points
#' @param method fmm, periodic, natural, hyman or monoH.FC (spline), scaledstineman, stineman or parabola (stinterp)
#' @return interpolated time series
#' @examples
#'   data("david_stoolA_otus")
#'   data("david_stoolA_metadata")
#'   days=david_stoolA_metadata[1,]
#'   sorted=sort(apply(david_stoolA_otus,1,sum),decreasing=TRUE,index.return=TRUE)
#'   davida.top=david_stoolA_otus[sorted$ix[1:100],]
#'   tsplot(interpolate(davida.top,time.vector=days))
#' @export

interpolate<-function(x, interval=NA, time.index=NA, time.vector=c(), method="stineman"){
  if(is.na(time.index) && length(time.vector)==0){
    stop("No time points given. Please provide either an index in the matrix or a vector with the time points.")
  }
  if(!is.na(time.index)){
    time.vector=x[time.index,]
  }
  time.vector=as.numeric(time.vector)
  # remove time points with missing values
  missing.indices=which(is.na(time.vector))
  if(length(missing.indices) > 0){
    good.samples=setdiff(c(1:ncol(x)),missing.indices)
    time.vector=time.vector[good.samples]
    x=x[,good.samples]
  }
  # no interval given: take the mode
  if(is.na(interval)){
    intervals=diff(time.vector)
    intervalcounts=table(intervals)
    interval=which(intervalcounts==max(intervalcounts))
    print(paste("Selected interval: ",interval,sep=""))
  }else{
    # make sure we have a positive integer as interval
    interval=round(abs(interval))
  }
  # time vector with regular intervals
  xout=seq(time.vector[1],time.vector[length(time.vector)],by=interval)
  print(paste("Length of time series: ",length(time.vector),sep=""))
  print(paste("Length of time series after interpolation: ",length(xout),sep=""))
  interpolated=matrix(NA,nrow=nrow(x),ncol=length(xout))
  # interpolate row by row
  for(i in 1:nrow(x)){
    if(method == "natural" || method == "fmm" || method == "periodic" || method == "hyman" || method == "monoH.FC"){
      interpolated[i,]=spline(time.vector,as.numeric(x[i,]),xout=xout, method=method)$y
    }else if(method == "scaledstineman" || method == "stineman" || method == "parabola"){
      interpolated[i,]=stinterp(time.vector,as.numeric(x[i,]),xout=xout, method=method)$y
    }else{
      stop("The given interpolation method is not supported.")
    }
  }
  # update time vector in the matrix if provided
  if(!is.na(time.index)){
    interpolated[time.index,]=xout
  }
  return(interpolated)
}
