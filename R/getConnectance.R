#' Compute connectance of an interaction matrix.
#'
#' The connectance is defined as E/Nˆ2-N, where E is the number
#' of realized arcs (the number of non-zero entries in the interaction matrix) 
#' and N*(N-1) the number of possible arcs. 
#' The diagonal (self-arcs) is excluded.
#'
#' @param A interaction matrix
#' @return the connectance
#' @examples
#' getConnectance(adjustc(generateA(N=10),c=0.5))
#' @export

getConnectance<-function(A){
	 N=nrow(A)
	 # exclude diagonal from observed and possible interactions
	 c=(length(A[A!=0])-N)/(ncol(A)*ncol(A)-N)
	 return(c)
 }