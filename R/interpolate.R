#' @title Time Series Interpolation
#' @description Wrapper function to interpolate a time series. If groups are specified, each
#' group is processed separately. The time points of each group are supposed to be in chronological order.
#'
#' @details If no interval is provided, the intervals in the time vector are computed and the most frequent one is taken as the interval.
#' The default interpolation method is stineman. Note that interpolation can introduce negative values in the abundance matrix.
#' @param x the time series matrix, rows are objects and columns are time points
#' @param interval the target intervals to be present after interpolation
#' @param time.index the row index holding time points
#' @param time.vector the vector holding time points
#' @param method fmm, periodic, natural, hyman or monoH.FC (spline), scaledstineman, stineman or parabola (stinterp)
#' @param groups a vector with group assignments and as many entries as there are samples
#' @param negzero set negative values (can be introduced by interpolation) to zero
#' @return interpolated time series
#' @examples
#'   data("david_stoolA_otus")
#'   data("david_stoolA_metadata")
#'   days=david_stoolA_metadata[1,]
#'   sorted=sort(apply(david_stoolA_otus,1,sum),decreasing=TRUE,index.return=TRUE)
#'   davida.top=david_stoolA_otus[sorted$ix[1:100],]
#'   tsplot(interpolate(davida.top,time.vector=days))
#' @export

interpolate<-function(x, interval=NA, time.index=NA, time.vector=c(), method="stineman", groups=c(), negzero=FALSE){
  if(is.na(time.index) && length(time.vector)==0){
    stop("No time points given. Please provide either an index in the matrix or a vector with the time points.")
  }
  if(length(groups)>0){
    if(length(groups)!=ncol(x)){
      stop("Each sample should have a group assigned.")
    }
  }
  if(!is.na(time.index)){
    time.vector=x[time.index,]
  }
  time.vector=as.numeric(time.vector)
  updated.time.vector=c()
  if(length(groups)>0){
    index=1
    unique.groups=unique(groups)
    interpolated=c()
    for(group in unique.groups){
      print(paste("Processing group",group))
      member.indices=which(groups==group)
      print(paste("Number of members",length(member.indices)))
      x.group=x[,member.indices]
      time.vector.sub=time.vector[member.indices]
      #print(time.vector.sub)
      temp=interpolateSub(x=x.group, time.index=time.index, time.vector=time.vector.sub, method=method, group=group, negzero=negzero)
      updated.time.vector=c(updated.time.vector,temp$time)
      if(index==1){
        interpolated=temp$interpolated
        index=2
      }else{
        interpolated=cbind(interpolated,temp$interpolated)
      }
    }
  }else{
    temp=interpolateSub(x=x, time.vector=time.vector, method=method)
    interpolated=temp$interpolated
    updated.time.vector=temp$time
  }
  #print(updated.time.vector)
  return(interpolated)
}

interpolateSub<-function(x, time.vector=c(), time.index=NA, interval=NA, method="stineman", group=NA, negzero=FALSE){
  # remove time points with missing values
  missing.indices=which(is.na(time.vector))
  if(length(missing.indices) > 0){
    print(paste("Removing",length(missing.indices),"samples with missing time point information."))
    good.samples=setdiff(c(1:ncol(x)),missing.indices)
    time.vector=time.vector[good.samples]
    x=x[,good.samples]
  }
  # no interval given: take the mode
  if(is.na(interval)){
    intervals=diff(time.vector)
    intervalcounts=table(intervals)
    print(intervalcounts)
    if(length(which(intervals<0))>0){
      stop("Time points must be provided in chronological order.")
    }
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
      interpolated.ts=spline(time.vector,as.numeric(x[i,]),xout=xout, method=method)$y
    }else if(method == "scaledstineman" || method == "stineman" || method == "parabola"){
      interpolated.ts=stinepack::stinterp(time.vector,as.numeric(x[i,]),xout=xout, method=method)$y
    }else{
      stop("The given interpolation method is not supported.")
    }
    if(negzero){
      neg.indices=which(interpolated.ts<0)
      interpolated.ts[neg.indices]=0
    }
    interpolated[i,]=interpolated.ts
  }
  # update time vector in the matrix if provided
  if(!is.na(time.index)){
    interpolated[time.index,]=xout
  }
  rownames(interpolated)=rownames(x)
  if(is.na(group)){
    colnames(interpolated)=xout
  }else{
    colnames.inter=c()
    for(tp in xout){
      colnames.inter=c(colnames.inter,paste(group,tp,sep="_"))
    }
    colnames(interpolated)=colnames.inter
  }
  res=list(interpolated,xout)
  # interpolated time series and interpolated time points
  names(res)=c("interpolated","time")
  return(res)
}
