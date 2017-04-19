#' @title Connectance
#' @description Compute connectance of an interaction matrix.
#' @details The connectance is defined as $E/N^2-N$, where $E$ is the number
#' of realized arcs (the number of non-zero entries in the interaction matrix)
#' and $N*(N-1)$ the number of possible arcs.
#' The diagonal (self-arcs) is excluded.
#'
#' @param A an interaction matrix
#' @return the connectance
#' @examples
#'   A <- cbind(c(-1,0,1),c(0,-1,0),c(-1,0,-1))
#'   x <- getConnectance(A)
#' @export
getConnectance <- function(A){

	 N <- nrow(A)
	 
	 # exclude diagonal from observed and possible interactions
	 c <- (length(A[A!=0])-N)/(ncol(A)*ncol(A)-N)
	 
	 return(c)
 }
