#' @title Rarefaction combined with sample filtering
#'
#' @description Rarefy a matrix to the given minimum count number column-wise
#' using vegan's rrarefy function. If columns have less than the minimum count number,
#' they are discarded.
#' @param x a matrix
#' @param min minimum count to which x is to be rarefied (if equal to zero, the minimum column sum is taken as min)
#' @return a list with the rarefied matrix (rar) and the indices of the columns that were kept (colindices)
#' @examples
#' data(david_stoolA_otus)
#' filtered=rarefyFilter(david_stoolA_otus,min=1000)
#' # the rarefied OTU table is stored in rar
#' stoolARar=filtered$rar
#' # print names of samples that were discarded because they had less than min reads
#' discarded=setdiff(1:ncol(david_stoolA_otus),filtered$colindices)
#' print(colnames(david_stoolA_otus)[discarded])
#' @export

rarefyFilter<-function(x,min = 0){
  keep=c()
  if(min < 0){
    stop("Min should be either 0 or positive.")
  }
  if(min == 0){
    min=min(colsums=apply(x,2,sum))
    print(paste("Rarefy to minimum count",min))
    keep=c(1:ncol(x))
  }else{
    colsums=apply(x,2,sum)
    # there are columns below the minimum
    if(min(colsums) < min){
      # loop column sums
      for(j in 1:ncol(x)){
        if(colsums[j] >= min){
          keep=c(keep,j)
        }
      }
      print(paste("Number of columns",ncol(x)))
      print(paste("Keeping ",length(keep)," columns with column sums equal or above",min))
      x=x[,keep]
    }
  }
  rar=t(vegan::rrarefy(t(x),min))
  res=list(rar,keep)
  names(res)=c("rar","colindices")
  return(res)
}
