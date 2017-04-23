#' @title Filter taxa in an abundance matrix
#' @description Discard taxa with less than the given minimum number of occurrences.
#' @param x taxon abundance matrix, rows are taxa, columns are samples
#' @param minocc minimum occurrence (minimum number of samples with non-zero taxon abundance)
#' @param keepSum If keepSum is true, the discarded rows are summed and the sum is added as a row with name: summed-nonfeat-rows
#' @return filtered abundance matrix
#' @examples
#' \dontrun{
#' data(david_stoolA_otus)
#' minocc=round(ncol(david_stoolA_otus)/3)
#' stoolAFiltered=filterTaxonMatrix(david_stoolA_otus,minocc=minocc)
#' print(paste("Filtered taxa with less than a minimum occurrence of:",minocc))
#' print(paste("Taxon number before filtering:",nrow(david_stoolA_otus))
#' print(paste("Taxon number after filtering:",nrow(stoolAFiltered)))
#' }
#' @export

filterTaxonMatrix<-function(x, minocc=0, keepSum=FALSE){
  toFilter=c()
  for(i in 1:nrow(x)){
    occurrences=0
    # count occurrences in row i
    for(j in 1:ncol(x)){
      if(!is.na(x[i,j]) && x[i,j]!=0){
        occurrences=occurrences+1
      }
    }
    if(occurrences<minocc){
      toFilter=append(toFilter,i)
    }
  }
  indices.tokeep=setdiff(c(1:nrow(x)),toFilter)
  #print(paste("Filtering",rownames(x)[toFilter]))
  if(keepSum==TRUE){
    filtered=x[toFilter,]
    x=x[indices.tokeep,]
    rownames=rownames(x)
    sums.filtered=apply(filtered,2,sum)
    x=rbind(x,sums.filtered)
    rownames=append(rownames,"summed-nonfeat-rows")
    rownames(x)=rownames
  }else{
    x=x[indices.tokeep,]
  }
  return(x)
}

