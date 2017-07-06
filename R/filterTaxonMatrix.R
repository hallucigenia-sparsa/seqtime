#' @title Filter taxa in an abundance matrix
#' @description Discard taxa with less than the given minimum number of occurrences.
#' @param x taxon abundance matrix, rows are taxa, columns are samples
#' @param minocc minimum occurrence (minimum number of samples with non-zero taxon abundance)
#' @param dependency if true, remove all taxa with a slope above -0.5 or a non-linear slope in the periodogram in log-scale (samples are supposed to represent equidistant time points)
#' @param keepSum If keepSum is true, the discarded rows are summed and the sum is added as a row with name: summed-nonfeat-rows
#' @param return.filtered.indices if true, return an object with the filtered abundance matrix in mat and the indices of removed taxa in the original matrix in filtered.indices
#' @return filtered abundance matrix
#' @examples
#' \dontrun{
#' data(david_stoolA_otus)
#' minocc=round(ncol(david_stoolA_otus)/3)
#' stoolAFiltered=filterTaxonMatrix(david_stoolA_otus,minocc=minocc)
#' print(paste("Filtered taxa with less than a minimum occurrence of:",minocc))
#' print(paste("Taxon number before filtering:",nrow(david_stoolA_otus)))
#' print(paste("Taxon number after filtering:",nrow(stoolAFiltered)))
#' }
#' @export

filterTaxonMatrix<-function(x, minocc=0, dependency=FALSE, keepSum=FALSE, return.filtered.indices=FALSE){
  toFilter=c()
  xcopy=x
  # convert into presence/absence matrix
  xcopy[xcopy>0]=1
  # sum for each taxon = number of occurrences across samples
  rowsums=apply(xcopy,1,sum)
  toFilter=which(rowsums<minocc)

  if(dependency==TRUE){
    # allow large deviations from pink and brown noise
    # minimum slope: -0.5
    nt=identifyNoisetypes(x, epsilon=0.5)
    # remove all non-pink, non-brown and non-black taxa
    toKeep=c(nt$pink, nt$brown, nt$black)
    toFilter=c(toFilter, setdiff(c(1:nrow(x)),toKeep))
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
  if(return.filtered.indices==TRUE){
    res=list(x,toFilter)
    names(res)=c("mat","filtered.indices")
    return(res)
  }else{
    return(x)
  }
}

