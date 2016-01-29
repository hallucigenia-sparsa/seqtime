#' Generate time series from Ricker model
#'
#' If an explosion occurs, -1 is returned. An explosion
#' is reported if any species crosses the explosion bound.
#'
#' @param N number of species
#' @param A interaction matrix
#' @param K carrying capacity
#' @param y initial abundances
#' @param sigma noise level
#' @param tend number of time points
#' @param tskip number of initial time points to skip (to skip the transient)
#' @param explosion.bound boundary for explosion
#' @return a matrix with species abundances as rows and time points as columns
#' @examples
#' tsplot(ricker(10,generateA(10),K=rep(0.01,10)),type="l", header="ricker")
#' @seealso \code{\link{simuntb}} for the neutral model and \code{\link{glv}} for the generalized Lotka Volterra model
#' @export

ricker<-function(N, A, K=rep(0.1,N), y=runif(N), sigma=0.01, tend=100, tskip=0, explosion.bound=10^4){
  out=matrix(nrow=N, ncol=tend-tskip)
  out[,1]=y
  # simulate difference equation
  for(t in 2:tend){
    b=rlnorm(N,meanlog=0,sdlog=sigma)
    y=b*y*exp(K*(1+A%*%y))
    if(max(y) > explosion.bound){
      print("Explosion!")
      return(-1)
    }
    else if(length(y[y<0]) > 0){
      stop("Species below 0!")
    }
    if(t > tskip){
      out[,t]=y
    }
  }
  return(out)
}

