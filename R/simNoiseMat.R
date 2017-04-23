#' @title Simulate Noise
#' @description Simulate a matrix composed of rows with different noise types.
#'
#' @details White noise is generated with rnorm and brown noise with cumsum(rnorm).
#' Pink noise is generated with the noise function from the tuneR package.
#' Brown noise can optionally be generated with geometric Brown motion via the
#' GBM function in the sde package. Option stayDead allows simulating extinction
#' events; i.e. once a value of zero or below has been hit, the row only takes zero values.
#'
#' @param samples the number of samples to be drawn
#' @param noisetypes a list specifying the number of rows to be generated for each noise type
#' @param brown.type the algorithm to generate brown noise, options: bm (Brown motion) and gbm (geometrical Brown motion; requires sde package)
#' @param mean parameter for white and brown noise (type bm)
#' @param sd parameter for white and brown noise (type bm)
#' @param c constant added to pink noise
#' @param r parameter for brown noise (type gbm)
#' @param stayDead once a time series hits zero or below, it stays at zero
#' @return a matrix where each row represents a distribution following a specified noise type
#' @examples
#'   # plot the power spectrum of pink noise
#'   ps <- powerspec(simNoiseMat(samples=100,noisetypes=list(white=0,brown=0,pink=1))[1,], plot=TRUE)
#'   # plot a matrix with mixed noise types
#'   mat <- simNoiseMat(samples=300,noisetypes=list(white=0,brown=20,pink=30))
#'   tsplot(mat,type="l", header="simulated noise")
#'   # Taylor law of brown noise
#'   t <- taylor(simNoiseMat(samples=300,noisetypes=list(white=0,brown=50,pink=0)),type="taylor")
#'   # Check generation of pink noise
#'   i <- identifyNoisetypes(simNoiseMat(samples=500,noisetypes=list(white=0,brown=0,pink=10)))
#' @export
simNoiseMat<-function(samples=100, noisetypes=list(white=2,brown=3,pink=5), brown.type="bm", mean=5, sd=1, c=1, r=2, stayDead=FALSE){
  noiseMat=matrix(nrow=noisetypes$white+noisetypes$brown+noisetypes$pink,ncol=samples)
  counter=1
  # generate requested number of instances of white noise
  if(noisetypes$white > 0){
    for(i in 1:noisetypes$white){
      # white noise:  slope in log-log spectral plot around 0
      noiseMat[counter,]=rnorm(samples,mean=mean,sd=sd);
      counter = counter+1
    }
  }
  # generate requested number of instances of brown noise
  if(noisetypes$brown > 0){
    # check whether sde is there
    if(brown.type == "gbm"){
      # https://stat.ethz.ch/pipermail/r-help/2005-September/078958.html
      searchSde=length(grep(paste("^package:","sde", "$", sep=""), search()))
      if (searchSde==0) {
          stop(paste("sde is not installed/loaded. It is required for option ",brown.type,".",sep=""))
      }
    }
    for(i in 1:noisetypes$brown){
      if(brown.type == "bm"){
        # brown noise:  slope in log-log spectral plot around -2
        noiseMat[counter,]=cumsum(rnorm(samples,mean=mean,sd=sd));
      }else if(brown.type == "gbm"){
        noiseMat[counter,]=sde::GBM(N=(samples-1),r=r)
      }else{
        stop("Please specify the algorithm for brown noise, either bm or gbm")
      }
      counter = counter+1
    }
  }
  # generate requested number of instances of pink noise
  if(noisetypes$pink > 0){
    for(i in 1:noisetypes$pink){
      pink=tuneR::noise(kind="pink")
      noiseMat[counter,]=pink@left[1:samples]
      counter = counter+1
    }
  }

  if(stayDead == TRUE){
    for(i in 1:nrow(noiseMat)){
      zeroThere = FALSE
      for(j in 1:ncol(noiseMat)){
        if(zeroThere == TRUE){
            noiseMat[i,j]=0
        }else{
          if(noiseMat[i,j] <=0){
            noiseMat[i,j]=0
            zeroThere = TRUE
          }
        }
      } # end column loop
    } # end row loop
  } # end stayDead

  return(noiseMat)
}

