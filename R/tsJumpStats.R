#' @title Compute statistics on jumps through community space
#'
#' @description Consecutive samples representing time points in an ordination plot can be interpreted
#' as vectors, which have a length and an angle. If ordinate is true, these lengths and angles
#' are computed. If ordinate is false, dissimilarities between consecutive time points are computed.
#' If a perturbation object is provided and plot.type hist is true, the histogram of jump lengths is
#' plotted and the significance of the difference between perturbed and non-perturbed jump lengths is
#' assessed with a Wilcoxon test.
#'
#' @param x a community time series matrix, with rows as taxa and columns as time points
#' @param time.given sample names provide time steps, only needed for plotting and plot.type jumps
#' @param time.unit unit of time, only needed for plotting and plot.type jumps
#' @param ordinate if TRUE, compute jumps in ordination plot
#' @param distance the distance (or dissimilarity) used to compute jumps directly or to compute the ordination
#' @param dimensions if ordinate is TRUE, the principal components considered
#' @param plot if TRUE, a plot of the jumps is made
#' @param plot.type jumps: plot jumps as bars in the order of occurrence, hist: a histogram of jump lengths for perturbed and non-perturbed sample pairs (only if perturbation object is provided)
#' @param header text to be used as plot title instead of default text
#' @param subsample for plot.type hist: subsample larger jump vector randomly down to smaller one, so that there are as many perturbed as non-perturbed jumps
#' @param perturb a perturbation object, if provided, the plot is colored accordingly
#' @examples
#' \dontrun{
#' data(david_stoolA_otus)
#' data(david_stoolA_metadata)
#' rarefaction=rarefyFilter(david_stoolA_otus,min = 10000)
#' rarefiedA=rarefaction$rar
#' days=david_stoolA_metadata[1,rarefaction$colindices]
#' interpA=interpolate(rarefiedA,time.vector=days,interval=1,method="stineman")
#' interpA[interpA<0]=0
#' perturbA=perturbation(times=c(80,105), durations=c(5,9))
#' res=tsJumpStats(interpA, plot=TRUE, perturb=perturbA, header="in stool A")
#' out=powerspec(res$lengths,plot=TRUE)
#' }
#' @export

tsJumpStats<-function(x, time.given=FALSE, time.unit="days", ordinate=TRUE, distance="bray", dimensions=c(1,2), plot=FALSE, plot.type="jumps", header="", subsample=FALSE, perturb=NULL){
  if(ordinate==TRUE){
    ordinate.res=vegan::capscale(data.frame(t(x))~1,distance=distance)
  }
  jumps=c()
  angles=c()
  for(sample.index in 1:(ncol(x)-1)){
    # vector defined by the two consecutive points in multidimensional space
    betweenvector=c()
    # vector defined by the null point and the first point
    firstvector=c()
    # vector defined by the null point and the second point
    secondvector=c()
    # compute the euclidean distance between points in multi-dimensional space
    if(ordinate==TRUE){
      # loop over dimensions
      for(dim.index in 1:length(dimensions)){
        firstvector=c(firstvector,ordinate.res$CA$u[sample.index,dimensions[dim.index]])
        secondvector=c(secondvector,ordinate.res$CA$u[(sample.index+1),dimensions[dim.index]])
        pointdim=ordinate.res$CA$u[(sample.index+1),dimensions[dim.index]] - ordinate.res$CA$u[sample.index,dimensions[dim.index]]
        betweenvector=c(betweenvector,pointdim)
      }
      betweenvector=betweenvector^2
      length=sqrt(sum(betweenvector))
      jumps=c(jumps,length)
      firstlength=sqrt(sum(firstvector^2))
      secondlength=sqrt(sum(secondvector^2))
      dotproduct=sum(firstvector*secondvector)
      angle=acos(dotproduct/(firstlength*secondlength))
      angles=c(angles,angle)
    }else{
      #print(sample.index)
      mat=rbind(x[,sample.index],x[,(sample.index+1)])
      jumps=c(jumps,vegdist(mat, distance=distance)[1])
      #jumps=c(jumps,distmat[sample.index,(sample.index+1)])
    }
  }
  if(plot==TRUE){
    if(plot.type=="hist" && !is.null(perturb)){
      perturb.indicator=getPerturbedIndices(1:ncol(x),perturb)
      perturb.indices=which(perturb.indicator==TRUE)
      normal.indices=which(perturb.indicator==FALSE)
      jump.perturb=jumps[perturb.indices]
      jump.normal=jumps[normal.indices]
      if(subsample==TRUE){
        if(length(jump.normal)>length(jump.perturb)){
            jump.normal=sample(jump.normal)[1:length(jump.perturb)]
        }else if(length(jump.perturb)>length(jump.normal)){
            jump.perturb=sample(jump.perturb)[1:length(jump.normal)]
        }
      }
      print(paste("Number of non-perturbed jumps:",length(jump.normal)))
      print(paste("Standard deviation of non-perturbed jumps:",sd(jump.normal,na.rm=TRUE)))
      print(paste("Number of perturbed jumps:",length(jump.perturb)))
      print(paste("Standard deviation of perturbed jumps:",sd(jump.perturb,na.rm=TRUE)))
      wilcox.out=wilcox.test(jump.normal,jump.perturb)
      print(wilcox.out)
      col2=rgb(0,1,0,0.5)
      col1=rgb(1,0,0,0.5)
      # limits
      xmax=max(jump.perturb, na.rm=TRUE)
      xmin=min(jump.perturb,na.rm=TRUE)
      ymax=max(jump.normal,na.rm=TRUE)
      ymin=min(jump.normal,na.rm=TRUE)
      max=max(xmax,ymax)
      min=min(ymin,xmin)
      xmaxD=max(hist(jump.perturb,breaks="FD",plot=FALSE)$density)
      ymaxD=max(hist(jump.normal,breaks="FD",plot=FALSE)$density)
      maxD=max(xmaxD,ymaxD)
      maxD=maxD+0.5 # add a margin
      if(header==""){
        title="Histogram of perturbed and non-perturbed jump lengths"
      }else{
        title=header
      }
      hist(jump.perturb,breaks="FD",xlim=c(min,max), ylim=c(0,maxD), prob=TRUE,col=col1, border=col1,xlab="Jump lengths", main=title)
      hist(jump.normal,breaks="FD",prob=TRUE,col=col2, border=col2,add=TRUE)
      legend("topright",legend=c("Perturbed","Normal"), lty = rep(1,2), col = c(col1,col2), merge = TRUE, bg = "white", text.col="black")
    }else{
      if(time.given==TRUE){
        time=as.numeric(colnames(x))
      }else{
        time=1:ncol(x)
      }
      defaultColor="green"
      perturbColor="red"
      colors=rep(defaultColor,length(time))
      if(!is.null(perturb)){
        perturb.indicator=getPerturbedIndices(1:ncol(x),perturb)
        perturb.indices=which(perturb.indicator==TRUE)
        normal.indices=which(perturb.indicator==FALSE)
        colors[perturb.indices]=perturbColor
      }
      colors=colors[2:(length(colors))]
      #print(length(colors))
      time=time[2:(length(time))]
      #print(length(time))
      # above 100 time points labels become unreadable
      if(length(time)<100){
        names(jumps)=time
      }
      if(header==""){
        title="Jumps in community composition space"
      }else{
        title=header
      }
      par(las=2, cex=0.9)
      barplot(jumps, col=colors, xlab=paste("Time",time.unit,sep=" in "), ylab="Dissimilarity", main=title)
    }
  }
  res=list(jumps,angles)
  names(res)=c("lengths","angles")
  res
}

# given a time series and a perturbation object, return a vector with FALSE for non-perturbed
# and true for perturbed samples
getPerturbedIndices<-function(time,perturb){
  perturbCounter=1
  durationCounter=1
  perturbationOn=FALSE
  indicator.timeseries=c()
  for(timepoint in time){
    applied=applyPerturbation(perturb=perturb,t=timepoint, perturbCounter=perturbCounter, durationCounter=durationCounter, perturbationOn=perturbationOn, ori.growthrates = c(), abundances=c())
    durationCounter=applied$durationCounter
    perturbCounter=applied$perturbCounter
    perturbationOn=applied$perturbationOn
    if(perturbationOn==TRUE){
      indicator.timeseries=c(indicator.timeseries,TRUE)
    }else{
      indicator.timeseries=c(indicator.timeseries,FALSE)
    }
  }
  return(indicator.timeseries)
}
