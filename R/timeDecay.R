#' @title Do a time-decay plot of a taxon matrix
#'
#' @description Plot the time against the sample distance.
#'
#' @param x a taxon matrix with rows representing taxa and columns samples
#' @param time vector with time steps
#' @param dissim sample dissimilarity to use, should be supported by vegdist
#' @param normtaxa divide each taxon vector in x by its sum
#' @param logdissim take the logarithm of the dissimilarity before fitting a line
#' @param logtime take the logarithm of the time interval before fitting a line
#' @param header a string to be appended to the plot title (Time decay)
#' @param units a string to describe the units of the time points (days, weeks etc)
#' @param plot whether to do the plot
#' @return a list with the time intervals, dissimilarity values, intersection, slope, p-value, adjusted R2, dissimilarity measure used and log status of dissimilarity measure
#' @examples
#' data("david_stoolA_otus")
#' data=rarefyFilter(david_stoolA_otus,min=10000)[[1]]
#' out.decay=timeDecay(data[,1:50], header="Stool subject A")
#' @export

timeDecay<-function(x, time=c(1:ncol(x)), dissim="bray", normtaxa=FALSE, logdissim=FALSE, logtime=FALSE, header="", units="", plot=TRUE){
  if(length(time) != ncol(x)){
    stop("The time vector has not as many entries as x has columns!")
  }
  if(normtaxa == TRUE){
    rowsums = apply(x,1,sum)
    # remove taxa with only zeros from matrix, to avoid dividing by a zero
    zero.indices=which(rowsums==0)
    rowsums=rowsums[setdiff(1:nrow(x),zero.indices)]
    x=x[setdiff(1:nrow(x),zero.indices),]
    x=x/rowsums
  }
  # compute sample-wise dissimilarities
  dissimMat=as.matrix(vegdist(t(x),method=dissim))
  intervals = c()
  dissimValues = c()
  # plot dissimilarities against time
  for(index1 in 1:(length(time)-1)){
    for(index2 in (index1+1):length(time)){
      interval = time[index2]-time[index1]
      # symmetric
      dissimVal = dissimMat[index1,index2]
      intervals=c(intervals,interval)
      dissimValues=c(dissimValues,dissimVal)
    }
  }
  # reg.data=data.frame(intervals,dissimValues)
  if(logdissim==TRUE){
    dissimValues=log(dissimValues)
  }
  if(logtime==TRUE){
    intervals=log(intervals)
  }
  if(length(dissimValues) > 1 && length(which(is.na(dissimValues)))==0 && length(which(is.infinite(dissimValues)))==0){
    linreg = lm(formula = dissimValues~intervals)
    intersection = linreg$coefficients[1]
    slope=linreg$coefficients[2]
    sum=summary(linreg)
    adjR2=sum$adj.r.squared
    pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
    #print(paste("slope",slope))
    #print(paste("p-value",pval))
    #print(paste("Adjusted R2:",sum$adj.r.squared))
    xlab="Interval"
    if(units != ""){
      xlab=paste(xlab,"in",units)
    }
    if(logdissim == TRUE){
      ylab=paste("Log community dissimilarity (",dissim,")", sep="")
    }else{
      ylab=paste("Community dissimilarity (",dissim,")", sep="")
    }
    if(logtime == TRUE){
      xlab="Log(interval)"
    }
    if(plot == TRUE){
      plot(intervals,dissimValues, xlab=xlab, ylab=ylab, main=paste("Time decay",header,"\nP-value",round(pval,3),", R2.adj",round(sum$adj.r.squared,3),", Slope",round(slope,3)))
      abline(linreg,bty="n",col="red")
    }
  }else{
    intersection=NA
    slope=NA
    pval=NA
    adjR2=NA
  }
  res=list(intervals, dissimValues, intersection, slope, pval, adjR2, dissim,logdissim)
  names(res)=c("intervals", "dissimvals","intersection","slope","pval","adjR2","dissim", "logdissim")
  return(res)
}
