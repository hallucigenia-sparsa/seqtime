#' @title Normalize a matrix
#'
#' @description Normalize a matrix column-wise by dividing each entry by its corresponding column sum.
#'
#' @details Columns summing to zero are removed by default.
#'
#' @param x a matrix
#' @param removeZero remove columns summing to zero.
#' @return a normalized matrix
#' @export

normalize<-function(x, removeZero=TRUE){
  colsums = apply(x,2,sum)
  # remove columns with only zeros from matrix, to avoid dividing by a zero
  if(removeZero==TRUE){
    zero.col.indices=which(colsums==0)
    colsums=colsums[setdiff(1:ncol(x),zero.col.indices)]
    x=x[setdiff(1:ncol(x),zero.col.indices),]
  }
  for(i in 1:ncol(x)){
    x[,i]=x[,i]/colsums[i]
  }
  x
}

