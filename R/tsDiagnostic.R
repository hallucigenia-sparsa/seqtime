#' Diagnostic plots for the outcome of a community simulation.
#'
#' @param x community time series, rows are species and columns are time points
#' @param type plot time vs number of species (spec), number of individuals (ind) or the distribution of the species with given index together with its proportion in the metacommunity, if provided
#' @param spec.index row index of selected species (required for type distrib)
#' @param spec.meta proportion of selected species in the meta-community (optional)
#' @example
#' N=50
#' m.vector=generateAbundances(N,mode=5,probabs=TRUE)
#' outH=simHubbell(N=N,I=1000,m.vector=m.vector,tend=1000)
#' tsDiagnostic(outH,type="distrib",spec.index = 4,spec.meta = m.vector[4])
#' tsDiagnostic(outH,type="spec")
#' tsDiagnostic(outH,type="ind")
#' @export

tsDiagnostic<-function(x, type="ind", spec.index=1, spec.meta=NA){
  if(type=="ind"){
    colsums=apply(x,2,sum)
    plot(1:ncol(x),colsums, ylab="Number of individuals", xlab="Time points") # no trend?
  }else if(type=="spec"){
    # species number bias - species number decreases with time
    specnum=c()
    for(j in 1:ncol(x)){
      v=as.numeric(x[,j])
      specnum=c(specnum, length(v[v>0]))
    }
    plot(1:ncol(x),specnum, ylab="Number of species", xlab="Time points")
  }else if(type=="distrib"){
    x=normalize(x)
    hist(x[spec.index,],main=paste("Distribution species",spec.index), xlab="Abundance")
    if(!is.na(spec.meta)){
      abline(v=spec.meta, col="blue")
    }
  }else{
    stop("Available plot types are ind, spec and distrib.")
  }
}
