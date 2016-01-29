#' Generate an interaction matrix
#'
#' Generate an interaction matrix, either randomly from a uniform distribution or
#' using Klemm-Eguiluz algorithm to generate a modular and scale-free interaction matrix.
#' @param type random (sample a uniform distribution), klemm (generate a Klemm-Eguiluz matrix) or empty (zero everywhere, except for diagonal which is set to d)
#' @param pep desired positive edge percentage (only for klemm)
#' @param d diagonal values (should be negative)
#' @param c desired connectance (interaction probability)
#' @param ignore.c do not adjust connectance
#' @param negedge.symm set symmetric negative interactions (only for klemm)
#' @param clique.size modularity parameter (only for klemm)
#' @return the interaction matrix
#' TODO: generate interaction strengths from half-normal distribution
#' @examples
#' klemm=generateA(N=10,type="klemm",c=0.5)
#' @references Klemm & Eguiluz, Growing Scale-Free Networks with Small World Behavior \url{http://arxiv.org/pdf/cond-mat/0107607v1.pdf}
#' @export

generateA<-function(N=100, type="random", c=0.02, ignore.c=FALSE, d=-0.5, pep=50, negedge.symm=FALSE, clique.size=5){
  A=matrix(0,nrow=N,ncol=N)      # init species interaction matrix
  if(type=="random"){
    for (i in 1:N){
      for(j in 1:N){
        if(i==j){
          A[i,j]=d
        }else{
          A[i,j] = runif(1,min=-0.5,max=0.5)
        }
      }
    }
  }else if(type=="empty"){
    diag(A)=d
  }else if(type=="klemm"){
    g<-klemm.game(N,verb=FALSE, clique.size)
    A=get.adjacency(g)
    A=as.matrix(A)
    diag(A)=d
  }

  if(ignore.c==FALSE){
    print(paste("Adjusting connectance to",c))
    A=adjustc(A,c=c)
  }

  # for Klemm-Eguiluz: inroduce negative edges and set random interaction strengths
  if(type=="klemm"){
    if(pep < 100){
      A=modifyA(A=A, perc=(100-pep), symmetric=negedge.symm, mode="negpercent")
    }
    # excluding edges on the diagonal
    print(paste("Final arc number (excluding self-arcs)", length(A[A!=0])-N ))
    # excluding negative edges on the diagonal
    print(paste("Final negative arc number (excluding self-arcs)", length(A[A<0])-N ))

    # check PEP and number of asymmetric negative interactions
    # assuming diagonal values are negative
    pep = getPep(A)
    print(paste("PEP:",pep))

    # convert binary interaction strengths (-1/1) into continuous ones using uniform distribution
    # zero would remove the edge, so the minimum strength is small, but non-zero
    min.strength=0.00001
    for(i in 1:nrow(A)){
      for(j in 1:nrow(A)){
        # skip diagonal
        if(i != j){
          A[i,j]=A[i,j]*runif(1,min=min.strength,max=1)
        }
      }
    }
  }
  return(A)
}



