#' @title Stability test for interaction matrix
#'
#' @description Test the stability of an interaction matrix.
#' The first two methods (coyte and eigen) determine whether the community's
#' steady state is stable in the sense that the community returns to it after
#' a small perturbation. The last method tests for explosion in simulations
#' with Ricker.
#'
#' @details Coyte's stability criterium is fulfilled if: max(r_e,r_s) - s < 0,
#' where s is the average of the diagonal values (the intra-species competition),
#' r_e is the half-horizontal radius of the eigenvalue ellipse of A
#' (which equals its Jacobian) and r_s is the eigenvalue corresponding to the average row sum.
#' The more negative max(r_e,r_s) - s is, the quicker the community returns
#' to steady state.
#' This criterium assumes that realized inter-species interaction strengths are drawn
#' from a half normal distribution |N(0,sigma2)| (the absolute of the normal
#' distribution). Only non-zero, non-diagonal interactions are considered to
#' compute the mean interaction strength.
#' The eigen method implements the classical stability criterion, which tests whether all real
#' parts of the eigenvalues of the interaction (=Jacobian) matrix are smaller than zero.
#'
#'
#' @param A the interaction matrix to test
#' @param method the method to test stability (coyte, eigen, ricker)
#' @param K carrying capacities for ricker test
#' @param y initial abundances for ricker test
#' @param sigma noise term for ricker test
#' @param explosion.bound explosion boundary for ricker test
#' @param tend end of simulation time for ricker test
#' @return boolean false if unstable and true if stable
#' @references Coyte et al. (2015). The ecology of the microbiome: Networks, competition, and stability. Science 350:663-666.
#' @examples
#' A=generateA(N=20)
#' testStability(A)
#' @export

testStability<-function(A, method="eigen", K=rep(0.1,N), y=runif(N), sigma=0.01, explosion.bound=10^4, tend=100){

  if(method == "coyte"){
    S = nrow(A)
    s=mean(diag(A)) # usually the diagonal value is a constant
    C = getConnectance(A) # excludes diagonal

    res=getAStats(A)
    # convert into proportions
    P=res$nbinteractions
    if(res$nbam>0){
      Pa=res$nbam/P
    }else{
      Pa=0
    }
    if(res$nbcomp>0){
      Pc=res$nbcomp/P
    }else{
      Pc=0
    }
    if(res$nbexp>0){
      Pe=res$nbexp/P
    }else{
      Pe=0
    }
    if(res$nbmut>0){
      Pm=res$nbmut/P
    }else{
      Pm=0
    }
    if(res$nbcom>0){
      Pp=res$nbcom/P
    }else{
      Pp=0
    }

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
    out=ricker(N, A=A, K=K, y=y, sigma=sigma, tend=tend, tskip=0, explosion.bound=explosion.bound)
    if(out[[1]]==-1){
      stable = FALSE
    }else{
      stable = TRUE
    }
  }

  return(stable)
}
