#' @title Assign Taxon Levels
#' @description Assign OTU identifiers or higher-level taxon names to an interaction matrix
#'
#' @param A an interaction matrix
#' @param data (optional) a matrix with taxon abundances, rows are taxa and columns are samples (only needed if type is non-empty)
#' @param metadata (optional) a matrix with metadata, samples match samples in OTU table (only needed if type is non-empty)
#' @param lineages a matrix with lineages, for the format please see data david_stool_lineages
#' @param type if non-empty, process data exactly as in method generateTS, supported: stoola and stoolb
#' @param taxon.level the taxon level to be assigned as row and column names, supported: otu, species, genus, family, order, class, phylum
#' @param uniqueNames make names unique by appending a counter if needed
#' @param higherLevelNames if given level is not known, assign the highest level that is known
#'
#' @examples
#' data("david_stool_lineages")
#' # make a random OTU interaction matrix as an example
#' N=30
#' A=generateA(N=N,c=0.1)
#' randomOTUIndices=sample(1:nrow(david_stool_lineages))[1:N]
#' rownames(A)=david_stool_lineages[randomOTUIndices,1]
#' colnames(A)=rownames(A)
#' Aclass=assignTaxonLevelsToA(A,lineages=david_stool_lineages, taxon.level="class")
#' classnetwork=plotA(Aclass,method="network")
#' @export

assignTaxonLevelsToA <- function(A, data=NULL, metadata=NULL, lineages, type="", taxon.level="genus", uniqueNames=FALSE, higherLevelNames=TRUE){

  # constants, only needed if type is non-empty
  david.minsamplesum=10000
  interpolation.method="stineman"
  N=nrow(A)

  if(type=="stoolb"){
    # omit the last sample in David data stool B, because there is a huge sampling gap of 66 days
    data=data[,1:(ncol(data)-1)]
    metadata=metadata[,1:(ncol(metadata)-1)]
  }
  if(type=="stoola" || type=="stoolb"){
    rarefyRes=rarefyFilter(data,min=david.minsamplesum)
    data=rarefyRes$rar

    # discard days with read count below minsamplesum
    days=metadata[1,rarefyRes$colindices]
    # make data equidistant
    stool.interp=interpolate(data,time.vector=days,interval=1,method=interpolation.method)
    sorted=sort(apply(stool.interp,1,sum),decreasing=TRUE,index.return=TRUE)
    otus=rownames(data)[sorted$ix[1:N]]
  }else{
    otus=rownames(A)
  }
  # match OTUs to entries in lineage table
  taxa=c()
  counter=1
  for(otu in otus){
    # ensure exact match
    if(taxon.level=="otu"){
      taxon=otu
    }else{
      index=which(lineages[,1]==otu)
      if(length(index)>0){
        if(taxon.level=="species"){
          levelIndex=8
        }else if(taxon.level=="genus"){
          levelIndex=7
        }else if(taxon.level=="family"){
          levelIndex=6
        }else if(taxon.level=="order"){
          levelIndex=5
        }else if(taxon.level=="class"){
          levelIndex=4
        }else if(taxon.level=="phylum"){
          levelIndex=3
        }else{
          stop("Requested taxon level not supported")
        }
        taxon=as.character(lineages[index,levelIndex])
      }else{
        warning("Did not find lineage of taxon",otu)
        taxon="unknown"
      }
    }
    # assign taxon name on level that is not none
    if(higherLevelNames==TRUE && taxon=="none"){
      levelIndex=levelIndex+1
      for(higherLevelIndex in levelIndex:2){
        if(as.character(lineages[index,higherLevelIndex])!="none"){
            taxon=as.character(lineages[index,higherLevelIndex])
            break
        }
      }
    }
    if(uniqueNames==TRUE){
      alreadySeen=which(taxa==taxon)
      if(length(alreadySeen)>0){
        taxon=paste(taxon,counter,sep="_")
        counter=counter+1
      }
    }
    taxa=append(taxa,taxon)
  }
  rownames(A)=taxa
  colnames(A)=rownames(A)
  return(A)
}
