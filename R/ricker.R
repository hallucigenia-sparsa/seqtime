#' @title Generate time series with the Ricker model
#'
#' @description The Ricker model is a discrete version of the generalized Lotka-Volterra model and
#' is implemented here as proposed by Fisher and Mehta in PLoS ONE 2014.
#'
#' @details If an explosion occurs, -1 is returned. An explosion
#' is reported if any species crosses the explosion bound. To avoid explosions, decrease the number of
#' interactions in the interaction matrix and/or increase the number of negative interactions in the
#' interaction matrix.
#'
#' @param N number of species
#' @param A interaction matrix
#' @param K carrying capacity
#' @param y initial abundances
#' @param sigma noise level, if set to a non-positive value, no noise is added
#' @param tend number of time points
#' @param tskip number of initial time points to skip (to skip the transient)
#' @param explosion.bound boundary for explosion
#' @param perturb a perturbation object
#' @return a matrix with species abundances as rows and time points as columns
#' @references Fisher & Mehta (2014). Identifying Keystone Species in the Human Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression. PLoS One 9:e102451
#' @examples
#' tsplot(ricker(10,generateA(10),K=rep(0.01,10)),header="ricker")
#' per=perturbation(times=seq(10,100,10),durations=rep(1,10),numberchanges=c(0.1,rep(0,9)))
#' tsplot(ricker(10,generateA(10),K=rep(0.01,10),perturb=per))
#' @seealso \code{\link{glv}} for the generalized Lotka Volterra model
#' @export
#'
#'

ricker<-function(N, A, K=rep(0.1,N), y=runif(N), sigma=0.05, tend=100, tskip=0, explosion.bound=10^8, perturb=NULL){
  out=matrix(nrow=N, ncol=tend-tskip)
  out[,1]=y
  perturbCounter=1
  durationCounter=1
  K.copy=K
  perturbationOn=FALSE
  # simulate difference equation
  for(t in 2:tend){
    if(sigma > 0){
      b=rlnorm(N,meanlog=0,sdlog=sigma)
    }else{
      b=rep(1,N)
    }
    if(!is.null(perturb)){
      applied=applyPerturbation(perturb,t=t,perturbCounter=perturbCounter,durationCounter=durationCounter,perturbationOn=perturbationOn,ori.growthrates=K.copy,abundances=y)
      y=applied$abundances
      K=applied$growthrates
      durationCounter=applied$durationCounter
      perturbCounter=applied$perturbCounter
      perturbationOn=applied$perturbationOn
      #print(perturbCounter)
      #print(perturbationOn)
    }
    y=b*y*exp(A%*%(y-K))
    if(max(y) > explosion.bound){
      # report which species explodes
      print("Explosion!")
      res=c(-1,which(y==max(y)))
      return(res)
    }
    else if(length(y[y<0]) > 0){
      stop("Species below 0!")
    }
    if(t > tskip){
      out[,t-tskip]=y
    }
  }
  return(out)
}



