#' @title Hill numbers
#'
#' @description The Hill numbers quantify biodiversity. The importance of the abundance distribution increases with increasing Hill order.
#' For q=0, the Hill number is the richness, for q=1, it is the (exponential) Shannon entropy and for q=2, it is the inverse Simpson index.
#' Note that the Hill order can also be a fraction, e.g. 0.5.
#'
#' @details The Hill number is defined as D=(SUM p_i^q)^1/(1-q), for i from 1 to S, where S is the species number,
#' p_i is the proportion of species i and q is the Hill order.
#' Since the Hill number involves a division by zero for q=1, please choose a sufficiently close q, such as 0.99999, when computing the
#' Hill number for 1.
#'
#' @param p relative abundance vector, should sum to one
#' @param q the Hill order
#' @return the Hill number
#' @references Hill (1973) "Diversity and evenness: a unifying notation and its consequences", Ecology 54: 427-432.
#' @examples
#' # even species distribution
#' hill(generateAbundances(N=1000,mode=6,k=0.001,probabs=TRUE),q=2)
#' # uneven species distribution
#' hill(generateAbundances(N=1000,mode=6,k=0.2,probabs=TRUE),q=2)
#' @export

hill<-function(p=rep(1/10,10), q=0){
  if(sum(p)!=1){
    warning("The relative abundances do not sum to one!")
  }
  if(q==1){
    stop("The Hill number involves a division by zero and is therefore not defined for q=1. Please choose a close q, e.g. 0.9999.")
  }
  if(q==0){
    S=length(which(p>0))
    return(S)
  }else{
    D=0
    for(i in 1:length(p)){
      D=D+p[i]^q
    }
    return(D^(1/(1-q)))
  }
}
