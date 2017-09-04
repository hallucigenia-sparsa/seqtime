#' @title Diagnostics for Community Simulation
#' @description Diagnostics for the outcome of a community simulation.
#' Reports the number of species that went extinct at the end of the simulation
#' and the number and proportion of zeros in the time series. In addition,
#' it offers diagnostic plots.
#' @param x community time series, rows are species and columns are time points
#' @param mode plot time vs number of species (spec), number of individuals/sum of abundance (ind) or the distribution of the species with given index (distrib)
#' @param plot plot the diagnostic plot of selected mode
#' @param spec.index row index of selected species (optional, required for mode distrib)
#' @param spec.meta proportion of selected species in the meta-community (optional, if provided it is added in plot mode distrib)
#' @param threshold threshold below which species are considered dead
#' @param K carrying capacity, in case of Hubbell and SOI, the number of individuals (optional, if provided it is shown in plot mode ind)
#' @param \\dots Additional arguments passed to plot()
#' @return depending on the mode, the function returns the number of surviving species over time, the number of individuals over time or the species distribution at the given time point
#' @examples
#' \dontrun{
#' N=50
#' m.vector=generateAbundances(N,mode=5,probabs=TRUE)
#' outH=simHubbell(M=N,N=N,I=1000,m.vector=m.vector,tend=1000)
#' diagnostics=tsDiagnostic(outH,mode="distrib",spec.index = 4,spec.meta = m.vector[4])
#' diagnostics=tsDiagnostic(outH,mode="spec")
#' diagnostics=tsDiagnostic(outH,mode="ind")
#' }
#' @export

tsDiagnostic<-function(x, mode="ind", plot=TRUE, spec.index=1, spec.meta=NA, threshold=0.01, K=NA, ...){
  numZeros=length(x[x==0])
  numValues=nrow(x)*ncol(x)
  result=c()
  onePerc=numValues/100
  percZeros=numZeros/onePerc
  lastTP=x[,ncol(x)]
  deadSpecies=which(lastTP<threshold)
  print(paste("Number of zeros:",numZeros))
  print(paste("Percentage of zeros:",percZeros))
  print(paste("Number of dead species at the final time point:",length(deadSpecies)))
  print("Indices of dead species at the final time point:")
  print(deadSpecies)
  if(mode=="ind"){
    colsums=apply(x,2,sum)
    result=colsums
    ylim=c(0,max(colsums))
    if(!is.na(K) && K > max(colsums)){
      ylim=c(0,K)
    }
    if(plot){
      plot(1:ncol(x),colsums,ylim=ylim,ylab="Sum of abundances", xlab="Time points", ...)
      if(!is.na(K)){
        abline(h=K, col="red")
      }
    }
  }else if(mode=="spec"){
    # species number over time
    specnum=c()
    for(j in 1:ncol(x)){
      v=as.numeric(x[,j])
      specnum=c(specnum, length(v[v>threshold]))
    }
    result=specnum
    if(plot){
      plot(1:ncol(x),specnum, ylab="Number of species", xlab="Time points", ...)
    }
  }else if(mode=="distrib"){
    if(spec.index < 0 || is.na(spec.index)){
      stop("Please provide the index of the species of interest!")
    }
    if(spec.index > nrow(x)){
      stop("The index of the species of interest is larger than the row number in the community time series!")
    }
    x=normalize(x, removeZero=FALSE)
    result=x[spec.index,]
    if(plot){
      hist(x[spec.index,],main=paste("Distribution species",spec.index), xlab="Abundance", ...)
      if(!is.na(spec.meta)){
        abline(v=spec.meta, col="blue")
      }
    }
  }else{
    stop("Available plot modes are ind, spec and distrib.")
  }
  return(result)
}
