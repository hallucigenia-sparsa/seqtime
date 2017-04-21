#' @title Positive Edge Percentage
#' @description Compute the positive edge percentage of an interaction matrix.
#' @details Values on the diagonal are ignored. More precisely, the positive arc percentage is computed, e.g.
#' the percentage of positive entries in the interaction matrix.
#' @param A the interaction matrix
#' @return the positive edge percentage
#' @examples
#'   A <- cbind(c(-1,0,1),c(0,-1,0),c(-1,0,-1))
#'   x <- getPep(A)
#' @export
getPep <- function(A){
  # excluding edges on the diagonal
  diag(A)=0
  num.pos=length(A[A>0])
  N=nrow(A)
  num.all=num.pos+length(A[A<0])
  one.perc=num.all/100
  pos.perc=num.pos/one.perc
  return(pos.perc)
}
