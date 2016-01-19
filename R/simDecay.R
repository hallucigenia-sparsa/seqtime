#' Plot sample-wise taxon dissimilarities against the abundance difference of selected taxa
#'
#' @param x taxon-sample matrix
#' @param taxa vector of taxon row indices or row names present in x
#' @param logdissim take the logarithm of the dissimilarity before fitting a line
#'
#' @examples
#' data=rarefyFilter(david_stoolA_otus,min=10000)
#' # Faecalibacterium OTUs
#' taxa=c("OTU_72853","OTU_119271")
#' out=simDecay(data,taxa)
#'
simDecay<-function(x,taxa=c(1),dissim="bray",logdissim=FALSE){
  if(length(taxa)==0){
    stop("Please provide at least the index or the name of one taxon.")
  }
  differences = c()
  dissimValues = c()
  taxonStr=taxa[1]
  # assemble printable taxon string
  for(taxon.index in 2:length(taxa)){
    taxonStr=paste(taxonStr,taxa[taxon.index],sep=", ")
  }
  if(is.character(taxa)){
    indices=c()
    # grep row indices from row names
    for(taxon.index in 1:length(taxa)){
      index=grep(paste(taxa[taxon.index],"\\b",sep=""),rownames(x),perl=TRUE)
      if(length(index)==0){
        print(paste("Taxon",taxa[taxon.index],"not found as a row name."))
      }else{
        indices=c(indices,index)
      }
    } # end loop over taxon names
    if(length(indices)==0){
      stop("Could not match any of the taxon names to row names!")
    }
    taxa=indices
  }
  abundances = x[taxa,]
  # compute sample-wise dissimilarities without selected taxa
  keep=setdiff(c(1:nrow(x)),taxa)
  x=x[keep,]
  dissimMat=as.matrix(vegdist(t(x),method=dissim))
  # plot dissimilarities against abundance difference
  for(index1 in 1:(ncol(x)-1)){
    for(index2 in (index1+1):ncol(x)){
      difference = sum(abs(abundances[,index2]-abundances[,index1]))
      # symmetric
      dissimVal = dissimMat[index1,index2]
      differences=c(differences,difference)
      dissimValues=c(dissimValues,dissimVal)
    }
  }
  if(logdissim==TRUE){
    dissimValues=log(dissimValues)
  }
  linreg = lm(formula = dissimValues~differences)
  intersection = linreg$coefficients[1]
  slope=linreg$coefficients[2]
  sum=summary(linreg)
  pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
  #print(paste("slope",slope))
  #print(paste("p-value",pval))
  #print(paste("Adjusted R2:",sum$adj.r.squared))
  xlab="Difference"
  if(length(taxa) > 1){
    xlab=paste(xlab,"in taxa",taxonStr)
  }else{
    xlab=paste(xlab,"in taxon",taxonStr)
  }
  plot(differences,dissimValues, xlab=xlab, ylab=paste("Community dissimilarity (",dissim,")", sep=""),main=paste("Dissimilarity change",header,"\nP-value",round(pval,3),", R2.adj",round(sum$adj.r.squared,3),", Slope",round(slope,3)))
  abline(linreg,bty="n",col="red")
  res=list(differences, dissimValues, intersection, slope, pval, sum$adj.r.squared, dissim,logdissim)
  names(res)=c("differences", "dissimvals","intersection","slope","pval","adjR2","dissim", "logdissim")
  return(res)
}
