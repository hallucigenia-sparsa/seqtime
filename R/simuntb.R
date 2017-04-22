#' @title Run the Unified Neutral Theory of Biodiversity (UNTB) model
#'
#' @description This function just provides a wrapper around the untb function in the untb package.
#'
#' @param N number of species
#' @param y initial abundances
#' @param m migration probability between 0 and 1
#' @param tskip number of initial time points to be skipped when returning the result (to avoid the transient)
#' @param tend number of time points (i.e. the number of generations)
#' @return a matrix with species abundances as rows and time points as columns
#' @seealso \code{\link{ricker}} for the Ricker model and \code{\link{glv}} for the generalized Lotka Volterra model
#' @examples
#' N=50
#' y=generateAbundances(N, mode=5)
#' tsplot(simUntb(N,y=y,tend=1000))
#' @export

simUntb<-function(N, y=rep(1,N), m=0.02, tskip=0, tend=5000){

  if(tend < tskip){
    stop("The total number of time points is smaller than the number of time points to be skipped!")
  }

  # gens: generations
  # keep: keep the whole time series
  # outcome: timepoints x species
  ts <- untb(start = y, prob = m, gens = tend, keep = TRUE, meta = untb::as.count(1:N))

  # skip the transient
  if(tskip > 0){
    ts <- ts[(tskip+1):nrow(ts),]
  }

  spec <- species.table(ts)

  return(t(spec))
}

