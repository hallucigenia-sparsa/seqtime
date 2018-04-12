#' @title Generate time series with the Ricker model
#'
#' @description The Ricker model is a discrete version of the generalized Lotka-Volterra model and
#' is implemented here as proposed by Fisher and Mehta in PLoS ONE 2014.
#'
#' @details If an explosion occurs, -1 is returned. An explosion
#' is reported if any species crosses the explosion bound. To avoid explosions, decrease the number of
#' interactions in the interaction matrix and/or increase the number of negative interactions in the
#' interaction matrix. K.trend contains carrying capacity percentages (from the original carrying capacity) that are added
#' to each species' current carrying capacity in each time step. K.trend allows to simulate an environmental trend.
#' Negative carrying capacities are set to zero. If both a perturbation and K.trend are specified, the perturbation is carried out first.
#'
#' @param N number of species
#' @param A interaction matrix
#' @param K a vector of carrying capacities
#' @param K.trend a vector of positive or negative carrying capacity percentages specified between 0 and 1
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

ricker<-function(N, A, K=rep(0.1,N), y=runif(N), sigma=0.05, K.trend=NA, tend=100, tskip=0, explosion.bound=10^8, perturb=NULL){
  if(length(y) != N){
    stop("y needs to have N entries.")
  }
  if(nrow(A)!=N || ncol(A)!=N){
    stop("A needs to have N rows and N columns.")
  }
  if(length(K)!=N){
    stop("K needs to have N entries.")
  }
  if(length(K.trend)>1 && length(K.trend)!=N){
    stop("K.trend needs to have N entries.")
  }
  if(tskip>=tend){
    stop("There are as many or more time points to skip than time points specified.")
  }
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
      if(perturb$times[1]==1){
        stop("Please do not specify a perturbation at the first time point.")
      }
      applied=applyPerturbation(perturb,t=t,perturbCounter=perturbCounter,durationCounter=durationCounter,perturbationOn=perturbationOn,ori.growthrates=K.copy,abundances=y)
      y=applied$abundances
      K=applied$growthrates
      durationCounter=applied$durationCounter
      perturbCounter=applied$perturbCounter
      perturbationOn=applied$perturbationOn
      #print(perturbCounter)
      #print(perturbationOn)
    }
    if(length(K.trend)>1){
      # calculate percentages to be added
      K.onepercent=K.copy/100
      K.percent=K.trend*100*K.onepercent
      K=K+K.percent
      # set negative carrying capacities to zero
      negative.K.indices=which(K<0)
      K[negative.K.indices]=0
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



