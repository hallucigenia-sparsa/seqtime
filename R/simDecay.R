#' @title Plot community similarity decay against selected taxa or metadata
#'
#' @description All pair-wise community dissimilarities are plotted
#' against the corresponding abundance difference of selected taxa or metadata.
#'
#' @param x taxon-sample matrix with taxa as rows and samples as columns
#' @param taxa vector of taxon row indices or row names present in x
#' @param metadata vector with values in the same order as samples in the taxon matrix
#' @param metadataName name of the metadata vector (e.g. age)
#' @param samplegroups supposed to assign an integer to each sample, with integers ranging from 1 to group number. If provided, intra-group dissimilarity dots are colored by sample group
#' @param dissim dissimilarity measure to use
#' @param header a string to be appended to the plot title (Dissimilarity change)
#' @param normtaxa divide each taxon vector in x by its sum
#' @param logdissim take the logarithm of the dissimilarity before fitting a line
#' @param intragroupsOnly only take sample pairs within the same group into account for dissimilarity computation
#'
#' @examples
#' N=50
#' data("david_stoolA_otus")
#' rar=rarefyFilter(david_stoolA_otus,min=10000)
#' data=rar[[1]]
#' # Faecalibacterium OTUs
#' taxa=c("OTU_72853","OTU_119271")
#' out1=simDecay(data[,1:N],taxa)
#' data("david_stoolA_metadata")
#' days=david_stoolA_metadata[1,rar[[2]]] # only keep samples that made it through rarefaction
#' out2=simDecay(data[,1:N],metadata=days[1:N],metadataName=rownames(david_stoolA_metadata)[1])
#'
#' @export

simDecay<-function(x,taxa=c(),metadata=c(), metadataName="", samplegroups=rep(1,ncol(x)), dissim="bray", normtaxa=FALSE, logdissim=FALSE, intragroupsOnly=FALSE, header=""){
  if(length(taxa)==0 && length(metadata) == 0){
    stop("Please provide at least the index or the name of one taxon or alternatively provide a metadata vector.")
  }
  if(length(metadata) > 0 && length(taxa) > 0){
    stop("Please provide either taxa or a metadata vector.")
  }
  if(length(metadata) > 0 && length(metadata) != ncol(x)){
    stop("There should be as many metadata entries as samples in the taxon matrix.")
  }
  if(length(samplegroups) != ncol(x)){
    stop("There should be as many sample group vector entries as samples in the taxon matrix.")
  }
  if(normtaxa == TRUE){
    rowsums = apply(x,1,sum)
    # remove taxa with only zeros from matrix, to avoid dividing by a zero
    zero.indices=which(rowsums==0)
    rowsums=rowsums[setdiff(1:nrow(x),zero.indices)]
    x=x[setdiff(1:nrow(x),zero.indices),]
    x=x/rowsums
  }
  differences = c()
  dissimValues = c()
  colors=c()
  groupNum=length(unique(samplegroups))
  colorlookup=list()
  if(groupNum == 1){
    hues=c("black")
  }else{
    generator = seq(0,1,1/groupNum)
    hues = hsv(generator)
  }
  taxonStr = ""
  if(length(taxa) > 0){
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
  }else{
    abundances=as.numeric(metadata)
  }
  dissimMat=as.matrix(vegdist(t(x),method=dissim))
  # plot dissimilarities against abundance difference
  for(index1 in 1:(ncol(x)-1)){
    for(index2 in (index1+1):ncol(x)){
      proceed=TRUE
      if(intragroupsOnly==TRUE){
        if(samplegroups[index1]!=samplegroups[index2]){
          proceed=FALSE
        }
      }
      if(proceed==TRUE){
        if(length(taxa) > 0){
          difference = sum(abs(abundances[,index2]-abundances[,index1]))
        }else{
          difference = abs(abundances[index2]-abundances[index1])
        }
        # symmetric
        dissimVal = dissimMat[index1,index2]
        differences=c(differences,difference)
        dissimValues=c(dissimValues,dissimVal)
        # color
        if(samplegroups[index1] == samplegroups[index2]){
          group=samplegroups[index1]
          col=hues[group]
          colorlookup[[paste("group",group)]]=col
          colors=c(colors,col)
        }else{
          colors=c(colors,"grey")
          #colorlookup[["intergroup"]]="grey"
        }
      }
    } # inner loop
  } # outer loop
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
  if(logdissim == TRUE){
    ylab=paste("Log community dissimilarity (",dissim,") difference", sep="")
  }else{
    ylab=paste("Community dissimilarity (",dissim,") difference", sep="")
  }
  xlab="Difference"
  if(length(taxa)>0){
    if(length(taxa) > 1){
      xlab=paste(xlab,"in taxa",taxonStr)
    }else{
      xlab=paste(xlab,"in taxon",taxonStr)
    }
  }else{
    xlab=paste(xlab,"in ",metadataName)
  }
  plot(differences,dissimValues, xlab=xlab, ylab=ylab, main=paste("Dissimilarity change",header,"\nP-value",round(pval,3),", R2.adj",round(sum$adj.r.squared,3),", Slope",round(slope,3)),col=colors)
  abline(linreg,bty="n",col="red")
  #lines(seq(0,1,by=0.1),seq(0,1,by=0.1),col="gray")
  if(groupNum > 1){
    legendnames=c()
    legendcolors=c()
    for(groupname in names(colorlookup)){
      if(groupname != ""){
        legendnames=c(legendnames,groupname)
        legendcolors=c(legendcolors,colorlookup[[groupname]])
      }
    }
    legend("topright",legend=legendnames,cex=0.9, bg = "white", text.col=legendcolors)
  }
  res=list(differences, dissimValues, intersection, slope, pval, sum$adj.r.squared, dissim,logdissim)
  names(res)=c("differences", "dissimvals","intersection","slope","pval","adjR2","dissim", "logdissim")
  return(res)
}
