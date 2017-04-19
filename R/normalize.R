#' Normalize a matrix
#'
#' Normalize a matrix column-wise by dividing each entry by its corresponding column sum
#'
#' @param x a matrix
#' @return a normalized matrix
#' @export

normalize<-function(x){
  colsums = apply(x,2,sum)
  for(i in 1:ncol(x)){
    x[,i]=x[,i]/colsums[i]
  }
  x
}

# Discard rows with less than the given minimum number of occurrences.
# If keep.sum is true, the discarded rows are added as a summed row with
# name: summed-nonfeat-rows
#
# x is a matrix with taxa as rows and samples as columns
#
filter<-function(x, minocc=0, keep.sum=TRUE){
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
  print(paste("Filtering",rownames(x)[toFilter]))
  if(keep.sum==TRUE){
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
