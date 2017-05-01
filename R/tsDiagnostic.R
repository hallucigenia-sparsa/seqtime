#' @title Diagnostics for Community Simulation
#' @description Diagnostics for the outcome of a community simulation.
#' Reports the number of species that went extinct at the end of the simulation
#' and the number and proportion of zeros in the time series. In addition,
#' it offers diagnostic plots.
#' @param x community time series, rows are species and columns are time points
#' @param type plot time vs number of species (spec), number of individuals/sum of abundance (ind) or the distribution of the species with given index (distrib)
#' @param plot plot the diagnostic plot of selected type
#' @param spec.index row index of selected species (optional, required for type distrib)
#' @param spec.meta proportion of selected species in the meta-community (optional, if provided it is added in plot type distrib)
#' @examples
#' \dontrun{
#' N=50
#' m.vector=generateAbundances(N,mode=5,probabs=TRUE)
#' outH=simHubbell(M=N,N=N,I=1000,m.vector=m.vector,tend=1000)
#' tsDiagnostic(outH,type="distrib",spec.index = 4,spec.meta = m.vector[4])
#' tsDiagnostic(outH,type="spec")
#' tsDiagnostic(outH,type="ind")
#' }
#' @export

tsDiagnostic<-function(x, type="ind", plot=TRUE, spec.index=1, spec.meta=NA){
  numZeros=length(x[x==0])
  numValues=nrow(x)*ncol(x)
  onePerc=numValues/100
  percZeros=numZeros/onePerc
  lastTP=x[,ncol(x)]
  deadSpecies=which(lastTP==0)
  print(paste("Number of zeros:",numZeros))
  print(paste("Percentage of zeros:",percZeros))
  print(paste("Number of dead species at the final time point:",length(deadSpecies)))
  print("Indices of dead species at the final time point:")
  print(deadSpecies)
  if(plot==TRUE){
    if(type=="ind"){
      colsums=apply(x,2,sum)
      plot(1:ncol(x),colsums, ylab="Sum of abundance", xlab="Time points") # no trend?
    }else if(type=="spec"){
      # species number over time
      specnum=c()
      for(j in 1:ncol(x)){
        v=as.numeric(x[,j])
        specnum=c(specnum, length(v[v>0]))
      }
      plot(1:ncol(x),specnum, ylab="Number of species", xlab="Time points")
    }else if(type=="distrib"){
      if(spec.index < 0 || is.na(spec.index)){
        stop("Please provide the index of the species of interest!")
      }
      if(spec.index > nrow(x)){
        stop("The index of the species of interest is larger than the row number in the community time series!")
      }
      x=normalize(x, removeZero=FALSE)
      hist(x[spec.index,],main=paste("Distribution species",spec.index), xlab="Abundance")
      if(!is.na(spec.meta)){
        abline(v=spec.meta, col="blue")
      }
    }else{
      stop("Available plot types are ind, spec and distrib.")
    }
  }
}
