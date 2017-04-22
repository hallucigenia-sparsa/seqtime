#' @title The DOC method by Bashan et al.
#'
#' @description Given  matrix of relative taxon abundances, plot the dissimilarity-overlap curve (DOC) sample-wise.
#' The DOC method was developed by Bashan and colleagues.
#'
#' @param x rows are taxa, columns are samples
#' @param B bootstrap iterations, if set to 0 or below, no bootstraps are carried out
#' @param polygons draw a polygon instead of the lower and upper confidence line
#' @param rand randomize the data sample-wise before starting the DOC analysis
#' @param lower.conf lower limit of the confidence interval
#' @param upper.conf upper limit of the confidence interval
#' @param null.model the null model to use, permut shuffles x sample-wise, assembly selects for each non-zero taxon one of the values taken across the samples at random
#' @return a list with the overlaps, dissimilarities, lowess smoothed overlaps and dissimilarites and lower and upper confidence intervals
#' @references A. Bashan et al. (2016). Universality of human microbial dynamics, Nature 534, 259-262.
#' @examples
#' \dontrun{
#' data("david_stoolA_otus")
#' data=rarefyFilter(david_stoolA_otus,min=10000)[[1]]
#' out=doc(data[,1:100],B=10) # apply the DOC method on the first 100 time points
#' N=50
#' M=500
#' metapop=generateAbundances(N=M, mode=5, probabs=TRUE)
#' ts=simHubbell(N=N, M=M,I=1500,d=N, m.vector=metapop, tskip=500, tend=600)
#' out2=doc(ts,B=10) # apply the DOC method to the Hubbell time series
#' }
#' @export

doc<-function(x, B=100, polygons=FALSE, rand=FALSE, lower.conf=0.03, upper.conf=0.97, null.model="assembly"){
  if(rand == TRUE){
    if(null.model == "permut"){
      # permute each sample separately
      for(i in 1:ncol(x)){
        x[,i]=sample(x[,i], replace=FALSE)
      }
    }else if(null.model == "assembly"){
      # keep the assembly: replace non-zero taxa by a randomly selected value from their range
      for(taxa in 1:nrow(x)){
        taxonvec=x[taxa,]
        nonzeroindices=which(taxonvec>0)
        for(samples in 1:ncol(x)){
          if(x[taxa,samples] > 0){
            x[taxa,samples]=taxonvec[sample(nonzeroindices,1)]
          } # taxon has a non-zero count
        } # sample loop
      } # taxon loop
      # renormalize sample sums to 1
      colsums = apply(x,2,sum)
      for(i in 1:ncol(x)){
        x[,i]=x[,i]/colsums[i]
      }
    }else{
      stop("Please choose null model permut or assembly.")
    }
  } # randomize
  res=computeOandD(x)
  comparisons=0.5*ncol(x)*(ncol(x)-1)
  bootstrapsO=matrix(NA,nrow=B,ncol=comparisons)
  bootstrapsD=matrix(NA, nrow=B, ncol=comparisons)
  lowerConfO=c()
  upperConfO=c()
  lowerConfD=c()
  upperConfD=c()
  if(B > 0){
    # do B bootstraps
    for(iteration in 1:B){
      print(paste("Bootstrap",iteration))
      # sample with replacement
      selected=sample(1:ncol(x),replace=TRUE)
      bootX=x[,selected]
      resB=computeOandD(bootX)
      outLB=lowess(resB$overlap,resB$dissim)
      bootstrapsO[iteration,]=outLB$x
      bootstrapsD[iteration,]=outLB$y
    }
    # get confidence intervals
    for(i in 1:comparisons){
      lowerConfO=c(lowerConfO, quantile(bootstrapsO[,i],lower.conf))
      upperConfO=c(upperConfO, quantile(bootstrapsO[,i],upper.conf))
      lowerConfD=c(lowerConfD, quantile(bootstrapsD[,i],lower.conf))
      upperConfD=c(upperConfD, quantile(bootstrapsD[,i],upper.conf))
    }
  }
  outL=lowess(res$overlap,res$dissim)
  if(polygons == FALSE){
    plot(res$overlap, res$dissim, type="p", pch=".", xlab="overlap", ylab="dissimilarity (rJSD)", main="DOC", col="darkgrey")
    # Lowess line
    lines(outL, type="l", col="blue")
    if(B>0){
      # lower confidence interval
      lines(lowerConfO, lowerConfD, type="l", col="red")
      # upper confidence interval
      lines(upperConfO, upperConfD, type="l", col="red")
    }
  }else{
    lowesscol="blue"
    fillcol=rgb(0,0,1,alpha=0.5)
    if(rand == TRUE){
      fillcol=rgb(1,0,0,alpha=0.5)
      lowesscol="red"
    }
    plot(res$overlap, res$dissim, type="n", pch=".", xlab="overlap", ylab="dissimilarity (rJSD)", main="DOC", col="darkgrey")
    # Lowess line
    lines(outL, type="l", col=lowesscol, lwd=2)
    if(B>0){
      polygon(upperConfO,upperConfD, col=fillcol, border="white")
      polygon(lowerConfO, lowerConfD, col="white", border="white")
    }
  }
  res=list(res$overlap, res$dissim, outL$x, outL$y, lowerConfO, lowerConfD, upperConfO, upperConfD)
  names(res)=c("overlap","dissim","lowesso","lowessd", "lowconfo", "lowconfd", "highconfo", "highconfd")
  return(res)
}

# rJSD of a single shared species is 0, since renormalization
# will lead to this species having an abundance of 1 in both samples.
# Thus, minShared is set to 2.
computeOandD<-function(x, minShared=2){
  overlaps=c()
  dissimilarities=c()
  # loop all sample pairs
  for(i in 1:(ncol(x)-1)){
    for(j in (i+1):ncol(x)){
      a=as.numeric(x[,i])
      b=as.numeric(x[,j])
      # shared species
      nonnulla=which(a>0)
      nonnullb=which(b>0)
      indices.shared=intersect(nonnulla, nonnullb)
      # sample pair shares at least the minimal number of taxa
      if(length(indices.shared) >= minShared){
        # overlap
        o=(sum(a[indices.shared])+sum(b[indices.shared]))/2
        # dissimilarity
        renorma=a[indices.shared]/sum(a[indices.shared])
        renormb=b[indices.shared]/sum(b[indices.shared])
        dissim=rJSD(renorma, renormb)
      # no taxon shared: overlap is zero, but rJSD is not defined (division by zero)
      # for relative abundances, rJSD is bounded by [0, sqrt(log(2))]
      # both are set to NA to skip the corresponding point
      }else{
        o = NA
        dissim = NA
      }
      overlaps=c(overlaps,o)
      dissimilarities=c(dissimilarities, dissim)
    }
  }
  res=list(overlaps, dissimilarities)
  names(res)=c("overlap","dissim")
  return(res)
}

# taken from: http://enterotype.embl.de/enterotypes.html
rJSD<-function(x,y){
 return(sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)))
}
# addition of a pseudocount, to avoid negative infinities
KLD<-function(x,y, pseudocount=10^-10){
  x = x + pseudocount
  y = y + pseudocount
  return(sum(x * log(x/y)))
}
