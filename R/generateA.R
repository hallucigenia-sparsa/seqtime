#' @title Generate an interaction matrix
#'
#' @description Generate an interaction matrix, either randomly from a uniform distribution or
#' using Klemm-Eguiluz algorithm to generate a modular and scale-free interaction matrix.
#'
#' @param N number of species
#' @param type random (sample a uniform distribution), klemm (generate a Klemm-Eguiluz matrix) or empty (zero everywhere, except for diagonal which is set to d)
#' @param pep desired positive edge percentage (only for klemm)
#' @param d diagonal values (should be negative)
#' @param min.strength random: minimal off-diagonal interaction strength (only for random)
#' @param max.strength random: maximal off-diagonal interaction strength, klemm: maximal absolute off-diagonal interaction strength
#' @param c desired connectance (interaction probability)
#' @param ignore.c do not adjust connectance
#' @param negedge.symm set symmetric negative interactions (only for klemm)
#' @param clique.size modularity parameter (only for klemm)
#' @param groups vector of group memberships for each species, assign NA if species does not belong to any group (only for random)
#' @param intra.group.strength interaction strength between members of the same group
#' @param inter.group.strength interaction strength between members of different groups (if not defined, will be assigned randomly)
#' @return the interaction matrix
#' @examples
#' klemm=generateA(N=10,type="klemm",c=0.5)
#' groups=c(rep(NA,5),rep(1,10),rep(2,5),rep(3,10),rep(4,10))
#' Agroup=generateA(N=40,groups=groups,c=0.5,intra.group.strength=0.1,inter.group.strength=-0.5, d=-1)
#' @references Klemm & Eguiluz, Growing Scale-Free Networks with Small World Behavior \url{http://arxiv.org/pdf/cond-mat/0107607v1.pdf}
#' @export

generateA<-function(N=100, type="random",pep=50, d=-0.5, min.strength=-0.5, max.strength=0.5, c=0.02, ignore.c=FALSE, negedge.symm=FALSE, clique.size=5, groups=c(), intra.group.strength=0.5, inter.group.strength=NA){
  A=matrix(0,nrow=N,ncol=N)      # init species interaction matrix
  if(type=="random"){
    if(length(groups)>0){
      if(length(groups)!=nrow(A)){
        stop("Please define a group membership for each species.")
      }
    }
    for (i in 1:N){
      for(j in 1:N){
        if(i==j){
          A[i,j]=d
        }else{
          if(length(groups)==0){
            A[i,j] = runif(1,min=min.strength,max=max.strength)
          }else{
            group1=groups[i]
            group2=groups[j]
            if(!is.na(group1) && !is.na(group2) && group1==group2){
              A[i,j] = intra.group.strength
            }else{
              # assign interaction strength between groups randomly
              if(is.na(inter.group.strength)){
                A[i,j] = runif(1,min=min.strength,max=max.strength)
              }else{
                # assign selected interaction strength between groups
                A[i,j] = inter.group.strength
              }
            }
          }
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
    A=modifyA(A,c=c, mode="adjustc")
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
    min.klemm.strength=0.00001
    for(i in 1:nrow(A)){
      for(j in 1:nrow(A)){
        # skip diagonal
        if(i != j){
          A[i,j]=A[i,j]*runif(1,min=min.klemm.strength,max=max.strength)
        }
      }
    }
  }
  return(A)
}



