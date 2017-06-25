#' @title Compute statistics on jumps through community space
#'
#' @description Consecutive samples representing time points in an ordination plot can be interpreted
#' as vectors, which have a length and an angle. If ordinate is true, these lengths and angles
#' are computed. If ordinate is false, dissimilarities between consecutive time points are computed.
#'
#' @param x a community time series matrix, with rows as taxa and columns as time points
#' @param ordinate if TRUE, compute jumps in ordination plot
#' @param distance if ordinate is FALSE, the distance used to compute jumps
#' @param ordinate.distance the distance used for ordination
#' @param dimensions if ordinate is TRUE, the principal components considered
#' @param plot if TRUE, a plot of the jumps is made
#' @param header text attached to plot title
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

tsJumpStats<-function(x, ordinate=TRUE, distance="bray", ordinate.distance="bray", dimensions=c(1,2), plot=FALSE, header="", perturb=NULL){
  if(ordinate==TRUE){
    ordinate.res=vegan::capscale(data.frame(t(x))~1,distance=ordinate.distance)
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
    defaultColor="black"
    perturbColor="red"
    colors=c(defaultColor)
    if(!is.null(perturb)){
      perturbCounter=1
      durationCounter=1
      perturbationOn=FALSE
      perturbIndices=c()
      for(timepoint in 1:ncol(x)){
        applied=applyPerturbation(perturb=perturb,t=timepoint, perturbCounter=perturbCounter, durationCounter=durationCounter, perturbationOn=perturbationOn, ori.growthrates = c(), abundances=c())
        durationCounter=applied$durationCounter
        perturbCounter=applied$perturbCounter
        perturbationOn=applied$perturbationOn
        if(perturbationOn==TRUE){
          colors=c(colors,perturbColor)
          perturbIndices=c(perturbIndices,timepoint)
        }else{
          colors=c(colors,defaultColor)
        }
      }
    }
    plot(jumps, type="b", pch="+",col=colors, xlab="Time", ylab="Dissimilarity", main=paste("Jumps in community configuration space",header))
  }
  res=list(jumps,angles)
  names(res)=c("lengths","angles")
  res
}
