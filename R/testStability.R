#' Stability test for interaction matrix
#'
#' Test the stability of an interaction matrix using Coyte's criterium
#' or the eigenvalue criterium.
#' Coyte's criterium is: max(r_e,r_s) - s < 0,
#' where s is the average of the diagonal values (the intra-species competition),
#' r_e is the half-horizontal radius of the eigenvalue ellipse of A and
#' r_s is the eigenvalue corresponding to the average row sum.
#' The more negative max(r_e,r_s) - s is, the quicker the community returns
#' to equilibrium.
#' This criterium assumes that realized inter-species interaction strengths are drawn
#' from a half normal distribution |N(0,sigma2)| (the absolute of the normal
#' distribution). Only non-zero, non-diagonal interactions are considered to
#' compute the mean interaction strength.
#' The eigenvalue criterium tests whether all real parts of the eigenvalues
#' of the interaction matrix are smaller than zero.
#' Ricker simply runs a simulation with Ricker and tests whether an explosion occurs.
#'
#'
#' @param A the interaction matrix to test
#' @param method the method to test stability (coyte, eigen, ricker)
#' @param K carrying capacities for ricker test
#' @param y initial abundances for ricker test
#' @param sigma noise term for ricker test
#' @param explosion.bound explosion boundary for ricker test
#' @return boolean false if unstable and true if stable
#' @references Coyte et al. (2015). The ecology of the microbiome: Networks, competition, and stability. Science 350:663-666.
#' @examples
#' A=generateA(N=20)
#' testStability(A)
#' @export

testStability<-function(A, method="eigen", K=rep(0.1,N), y=runif(N), sigma=0.01, explosion.bound=10^4){

  if(method == "coyte"){
    S = nrow(A)
    s=mean(diag(A)) # usually the diagonal value is a constant
    C = getConnectance(A) # excludes diagonal

    res=getAStats(A)
    # convert into proportions
    P=res$nbinteractions
    Pa=res$nbam/P
    Pc=res$nbcomp/P
    Pe=res$nbexp/P
    Pm=res$nbmut/P
    Pp=res$com/P

    var=res$varstrength
    EX=res$meanstrength
    EXsquare=EX*EX

    # compute stability criterion
    re=sqrt(S*C* (var*(1-((Pp+Pa)/2))-C*EXsquare*(Pm-Pc+((Pp-Pa)/2))^2) )
    re = re * (1+ EXsquare*( (2*Pm + 2*Pc + Pp + Pa -1) - C*(Pm - Pc + 0.5*Pp - 0.5*Pa)^2 )/ (var*(1-0.5*Pp - 0.5*Pa) - C*EXsquare*(Pm - Pc + 0.5*Pp - 0.5*Pa)^2) )

    rs = (S-1)*C*( Pm-Pc+((Pp-Pa)/2) )*EX

    stable = FALSE

    print(paste("re:",re))
    print(paste("rs:",rs))
    print(paste("s:",s))

    if(max(re,rs) - abs(s) < 0){
      stable = TRUE
    }
  }else if(method == "eigen"){
    eig=eigen(A)
    real=Re(eig$values)
    if(length(real[real>0]) > 0){
      stable = FALSE
    }else{
      stable = TRUE
    }
  }else if(method == "ricker"){
    N=nrow(A)
    out=ricker(N, A=A, K=K, y=y, sigma=sigma, tend=100, tskip=0, explosion.bound=explosion.bound)
    if(out[[1]]==-1){
      stable = FALSE
    }else{
      stable = TRUE
    }
  }

  return(stable)
}
