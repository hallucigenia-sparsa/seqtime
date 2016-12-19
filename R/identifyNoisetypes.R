#' Identify noise types in a matrix row-wise
#'
#' Row sums need to be above the given abundance threshold and the spectral density
#' needs to scale significantly with the frequency (p-value below 0.05) in log-log scale.
#' The periodogram is computed with stats::spectrum.
#' The function returns a noisetypes object, which groups matrix row indices
#' by noise type.
#' @param x a matrix with objects as rows and time points as columns
#' @param pval.threshold significance threshold for periodogram powerlaw goodness of fit
#' @param abund.threshold minimum sum per row
#' @return S3 noisetypes object
#' @examples
#' N=10
#' ricker.out=ricker(N,generateA(N),K=rep(0.01,N))
#' noisetypes=identifyNoisetypes(ricker.out)
#' plot(ricker.out[noisetypes$brown[1],], main=paste("Simulated OTU",noisetypes$brown[1]),ylab="Abundance")
#' @export

identifyNoisetypes<-function(x, pval.threshold = 0.05, abund.threshold=10){
  pink=c()
  brown=c()
  black=c()
  white=c()
  belowThreshold=c()
  nonclass=c()
  slopes.nonclass=c()
  nonsig=c()
  # loop over rows of x
  for(i in 1:nrow(x)){
    sum.taxon=sum(x[i,])
    # total abundance should be above threshold
    if(sum.taxon > abund.threshold){
      # spectrum is hidden by igraph, package name required
      out=stats::spectrum(x[i,], plot=FALSE)
      # avoid zeros
      if(length(which(out$freq==0))==0 && length(which(out$spec==0))==0){
        logfreq=log(out$freq)
        logspec=log(out$spec)
        reg.data=data.frame(logfreq,logspec)
        linreg = lm(formula = logspec~logfreq)
        slope=linreg$coefficients[2]
        sum=summary(linreg)
        pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
        if(pval < pval.threshold){
          # white noise: [-0.2, 0.2]
          if(slope < 0.2 && slope > -0.2){
            white=c(white,i)
            # pink noise:  [-1.2,-0.8]
          }else if(slope < -0.8 && slope > -1.2){
            pink=c(pink,i)
            # brown noise: [-2.2, -1.5]
          }else if(slope < -1.5 && slope > -2.2){
            brown=c(brown,i)
            # black noise
          }else if(slope < -3){
            black=c(black,i)
          }else{
            nonclass=c(nonclass,i)
            slopes.nonclass=c(slopes.nonclass,slope)
          }
        }else{
          nonsig=c(nonsig,i)
        }
      }else{
        nonsig=c(nonsig,i)
      }
    } # taxon too rare
    else{
      belowThreshold=c(belowThreshold,i)
    }
  } # end taxon loop
  print(paste("Number of taxa below the abundance threshold: ",length(belowThreshold)))
  print(paste("Number of taxa with non-significant power spectrum laws: ",length(nonsig)))
  print(paste("Number of taxa with non-classified power spectrum: ",length(nonclass)))
  print(paste("Number of taxa with white noise: ",length(white)))
  print(paste("Number of taxa with pink noise: ",length(pink)))
  print(paste("Number of taxa with brown noise: ",length(brown)))
  print(paste("Number of taxa with black noise: ",length(black)))
  noisetypes=list(white,pink,brown,black,nonclass,nonsig, slopes.nonclass,belowThreshold)
  names(noisetypes)=c("white","pink","brown","black","nonclass","nonsig","slopes.nonclass","belowT")
  class(noisetypes) <- "noisetypes"
  return(noisetypes)
}
