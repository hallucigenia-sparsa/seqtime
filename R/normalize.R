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
