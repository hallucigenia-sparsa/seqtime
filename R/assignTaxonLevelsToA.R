#' @title Assign Taxon Levels 
#' @description Assign OTU identifiers or higher-level taxon names to an interaction matrix inferred from the David data processed with method generateTS.
#'
#' @param A an interaction matrix inferred from the David data processed with method generateTS
#' @param type the data set, supported: stoola and stoolb
#' @param taxon.level the taxon level to be assigned as row and column names, supported: otu, species, genus, family, order, class, phylum
#' @param uniqueNames make names unique by appending a counter if needed
#' @param higherLevelNames if given level is not known, assign the highest level that is known
#'
#' @examples
#' stoolA <- assignTaxonLevelsToA(stoolA,type="stoola")
assignTaxonLevelsToA <- function(A, type="stoola", taxon.level="genus", uniqueNames=FALSE, higherLevelNames=TRUE){
  david.minsamplesum=10000
  interpolation.method="stineman"
  N=nrow(A)
  data("david_stool_lineages")
  if(type=="stoola"){
    data("david_stoolA_otus")
    data("david_stoolA_metadata")
    stool=david_stoolA_otus
    metadata=david_stoolA_metadata
  }else if(type=="stoolb"){
    data("david_stoolB_otus")
    data("david_stoolB_metadata")
    stool=david_stoolB_otus
    metadata=david_stoolB_metadata
  }
  if(type=="stoolb"){
    # omit the last sample in David data stool B, because there is a huge sampling gap of 66 days
    stool=stool[,1:(ncol(stool)-1)]
    metadata=metadata[,1:(ncol(metadata)-1)]
  }
  rarefyRes=rarefyFilter(stool,min=david.minsamplesum)
  stool=rarefyRes$rar
  # discard days with read count below minsamplesum
  days=metadata[1,rarefyRes$colindices]
  # make data equidistant
  stool.interp=interpolate(stool,time.vector=days,interval=1,method=interpolation.method)
  sorted=sort(apply(stool.interp,1,sum),decreasing=TRUE,index.return=TRUE)
  otus=rownames(stool)[sorted$ix[1:N]]
  # match OTUs to entries in lineage table
  taxa=c()
  counter=1
  for(otu in otus){
    # ensure exact match
    if(taxon.level=="otu"){
      taxon=otu
    }else{
      index=which(david_stool_lineages[,1]==otu)
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
        taxon=as.character(david_stool_lineages[index,levelIndex])
      }else{
        warning("Did not find lineage of taxon",otu)
        taxon="unknown"
      }
    }
    # assign taxon name on level that is not none
    if(higherLevelNames==TRUE && taxon=="none"){
      levelIndex=levelIndex+1
      for(higherLevelIndex in levelIndex:2){
        if(as.character(david_stool_lineages[index,higherLevelIndex])!="none"){
            taxon=as.character(david_stool_lineages[index,higherLevelIndex])
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
