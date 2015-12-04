#' Identify noise types in a count matrix row-wise
#'
#' Row count sums need to be above the given abundance threshold and
#' the spectrum needs to scale significantly with the frequency below p-value 0.05
#' in log-log scale.
#' The periodogram is computed with spectrum.
#' @param x a count matrix
#' @param pval.threshold significance threshold for periodogram powerlaw goodness of fit
#' @param abund.threshold minimum count sum per row
#' @return S3 noisetypes object
#' @examples
#' N=10
#' ricker.out=ricker(N,generateA(N),K=rep(0.01,N))
#' identify.noisetypes(ricker.out)
#' @export

identifyNoisetypes<-function(x, pval.threshold = 0.05, abund.threshold=10){
  pink=c()
  brown=c()
  white=c()
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
        }else{
          nonclass=c(nonclass,i)
          slopes.nonclass=c(slopes.nonclass,slope)
        }
      }else{
        nonsig=c(nonsig,i)
      }
    } # taxon too rare
  } # end taxon loop
  print(paste("Number of taxa with non-significant power spectrum laws: ",length(nonsig)))
  print(paste("Number of taxa with non-classified power spectrum: ",length(nonclass)))
  print(paste("Number of taxa with white noise: ",length(white)))
  print(paste("Number of taxa with pink noise: ",length(pink)))
  print(paste("Number of taxa with brown noise: ",length(brown)))
  noisetypes=list(white,pink,brown,nonclass,nonsig, slopes.nonclass)
  names(noisetypes)=c("white","pink","brown","nonclass","nonsig","slopes.nonclass")
  class(noisetypes) <- "noisetypes"
  return(noisetypes)
}
