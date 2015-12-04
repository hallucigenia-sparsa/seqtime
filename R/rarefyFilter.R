#' Rarefaction combined with sample filtering
#'
#' Rarefy a given count vector or matrix column-wise using vegan's rrarefy function.
#' If columns have less than the minimum count number, they are discarded.
#' @param x a matrix or vector
#' @param min minimum count to which x is to be rarefied (if equal to zero, minimum column sum is taken as min)
#' @return rarefied vector or matrix
#' @export

rarefyFilter<-function(x,min = 0){
  if(min < 0){
    stop("Min should be either 0 or positive.")
  }
  if(min == 0){
    min=min(colsums=apply(x,2,sum))
    print(paste("Rarefy to minimum count",min))
  }else{
    colsums=apply(x,2,sum)
    # there are columns below the minimum
    if(min(colsums) < min){
      keep=c()
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
  rar
}
