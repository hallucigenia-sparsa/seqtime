#' Compute the positive edge percentage of the interaction matrix
#'
#' Values on the diagonal are ignored.
#'
#' @param the interaction matrix
#' @return the positive edge percentage
#' @examples
#' A=cbind(c(-1,0,1),c(0,-1,0),c(-1,0,-1))
#' getPep(A)
#' @export

getPep<-function(A){
  # excluding edges on the diagonal
  diag(A)=0
  num.pos=length(A[A>0])
  N=nrow(A)
  num.all=num.pos+length(A[A<0])
  one.perc=num.all/100
  pos.perc=num.pos/one.perc
  return(pos.perc)
}
