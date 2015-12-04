#' Compute the positive edge percentage of the interaction matrix
#'
#' @param the interaction matrix (values on the diagonal should be negative)
#' @return the positive edge percentage
#' @export

getPep<-function(A){
  num.pos=length(A[A>0])
  # excluding negative edges on the diagonal
  N=nrow(A)
  num.all=length(A[A!=0])-N
  one.perc=num.all/100
  pos.perc=num.pos/one.perc
  return(pos.perc)
}
