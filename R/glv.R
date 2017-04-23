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
#' @return a matrix with species abundances as rows and time points as columns, column names give time points
#' @seealso \code{\link{ricker}} for the Ricker model
#' @examples
#' tsplot(glv(N=4,generateA(4)),header="gLV",type="l")
#' @export

glv<-function(N=4,A,b=runif(N),y=runif(N),tstart=0,tend=100,tstep=0.1){
  # parms as matrix
  parms=cbind(b,A)
  parms=cbind(rep(N,N),parms)
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
  dydt <- y*(b+A %*% y)
  list(dydt)
}
