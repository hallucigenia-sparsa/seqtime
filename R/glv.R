#' @title Simulate time series with the generalized Lotka-Volterra model
#'
#' @description Simulate a community time series using the generalized Lotka-Volterra model, defined as
#' \eqn{\frac{dx}{dt}=x(b+Ax)}, where x is the vector of species abundances, A is the interaction matrix
#' and b the vector of growth rates.
#'
#' @param N species number
#' @param A interaction matrix
#' @param b growth rates
#' @param y initial abundances
#' @param tstart initial time point
#' @param tend final time point
#' @param tstep time step
#' @param perturb perturbation object describing growth changes
#' @return a matrix with species abundances as rows and time points as columns, column names give time points
#' @seealso \code{\link{ricker}} for the Ricker model
#' @examples
#' tsplot(glv(N=4,generateA(4)),header="gLV")
#' @export

glv<-function(N=4,A,b=runif(N),y=runif(N),tstart=0,tend=100,tstep=0.1, perturb=NULL){
  # checks perturb object and includes in parms
  # parms as matrix
  parms=cbind(b,A)
  parms=cbind(rep(N,N),parms)
  if (!is.null(perturb)){
    if (length(perturb$growthchanges) != N){
      stop("Please provide as many growth changes as species!")
    }
    else {
      growthchanges = perturb$growthchanges
      parms = cbind(parms,growthchanges)
      if (length(perturb$times) != 0){
        if (length(perturb$times) > N){
          stop("Number of perturbations needs to be less than the number of species!")
        }
        else {
          count = length(perturb$times)
          times = vector(length = N, mode="numeric")
          times[1:(length(perturb$times))] = perturb$times
          durations = vector(length = N, mode="numeric")
          durations[1:(length(perturb$durations))] = perturb$durations
          parms = cbind(parms, times, durations)
          parms = cbind(parms, rep(count, N))
        }
      }
    }
  }
  times<-seq(tstart, tend, by=tstep)
  # run the simulation
  commtime<-lsoda(y, times, glvsolve, parms)
  time=commtime[,1]
  commtime=commtime[,2:ncol(commtime)]
  commtime=t(commtime)
  colnames(commtime)=time
  return(commtime)
}

# ==============================================================
# Equations (generalized Lotka-Volterra)
# ==============================================================

# matrix formulation of the ODE set
# t: current simulation time
# y: vector with current values of state variables (initial conditions)
# parms: parameter values
#
glvsolve<-function(t, y, parms){
  N=parms[1,1]  # species number
  b=parms[,2]   # vector of growth rates
  A=parms[,3:(N+2)] # interaction matrix
  if (length(parms[1,]) == (N+3)){
    c=parms[,(N+3)]
    b = b + c
  }
  else if (length(parms[1,]) > (N+3)){
    c=parms[,(N+3)]
    count = parms[1,(N+6)]
    times = parms[,(N+4)]
    durations = parms[,(N+5)]
    intervals = list()
    for (i in 1:count){
      tstart = times[i]
      tend = times[i] + durations[i]
      if ((t>tstart) & (t<tend)){
        b = b + c
        #print(b)
      }
    }
  }
  dydt <- y*(b+A %*% y)
  list(dydt)
}

