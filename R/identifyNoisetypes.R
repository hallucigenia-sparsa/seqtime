#' @title Identify Noise Types
#' @description Identify noise types in a matrix row-wise.
#' @details Row sums need to be above the given abundance threshold and the spectral density needs to scale significantly with the frequency (p-value below 0.05) in log-log scale. The periodogram is computed with \code{\link{spectrum}} from the \pkg{stats} package. The function returns a noisetypes object, which groups matrix row indices by noise type.
#' @param x a matrix with objects as rows and time points as columns
#' @param epsilon allowed deviation from the expected slope of 0 for white noise, -1 for pink noise and -2 for brown noise (all rows with a slope below -3 are classified as having black noise)
#' @param permut permute time points before computing noise types
#' @param pval.threshold significance threshold for periodogram powerlaw goodness of fit
#' @param abund.threshold minimum sum per row
#' @return S3 noisetypes object
#' @examples
#' N <- 10
#' ricker.out <- ricker(N,generateA(N),K=rep(0.01,N))
#' noise <- identifyNoisetypes(ricker.out, abund.threshold=0)
#' plot(noise)
#'
#' @export
identifyNoisetypes <- function(x, epsilon = 0.2, pval.threshold = 0.05, permut=FALSE, abund.threshold=0){
  if(epsilon < 0 || epsilon > 0.5){
    stop("Please select a value between 0 and 0.5 for epsilon.")
  }
  pink=c()
  brown=c()
  black=c()
  white=c()
  belowThreshold=c()
  nonclass=c()
  slopes.nonclass=c()
  nonsig=c()
  if(permut==TRUE){
    indices=sample(1:ncol(x))
    x=x[,indices]
  }
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
          # white noise:
          if(slope < epsilon && slope > -epsilon){
            white=c(white,i)
            # pink noise:
          }else if(slope < -(1-epsilon) && slope > -(1+epsilon)){
            pink=c(pink,i)
            # brown noise:
          }else if(slope < -(2-epsilon) && slope > -(2+epsilon)){
            brown=c(brown,i)
            # black noise:
          }else if(slope <= -3){
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
  res=noisetypes(white,pink,brown,black)
  return(res)
}
