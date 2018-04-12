#' @title Perturbation
#' @description A perturbation object which defines the number of times a perturbation is applied,
#' its duration for each time, the effects on the growth rates and/or the abundances. The perturbation
#' is assumed to have each time the same effect on the species' growth rates and/or abundances.
#' Note that perturbations should be separated by at least one time point.
#' If perturbations are repeatedly applied, the growth rate changes will take effect only in the first
#' time point of the duration and will then set back to the original values at the end of the duration.
#' The abundance changes will be applied at each time point at the duration.
#' Specifically, growth changes in Ricker will be applied to the carrying capacity vector. In the neutral
#' and SOI models, which do not implement growth rates, they will be applied to the metacommunity composition
#' (neutral model) and immigration vector (SOI model), such that immigration of species with increased "growth rate"
#' is more likely during perturbation. For the neutral model, capacityConstant is always set to TRUE.
#' If capacityConstant is false in the SOI model, but the total abundance is larger than the number of individuals,
#' randomly selected individuals will be removed until the total abundance equals the number of individuals.
#'
#' @param times a vector of time points at which perturbations take place
#' @param durations a vector of the same length as times that specifies the duration of each perturbation
#' @param growthchanges (optional) a vector of the same length as species in the (meta-)community that describes by what amount each species growth rate is decreased or increased (zero: no change, minus: decrease, else increase)
#' @param deathrate (optional) a constant (neutral model) or a vector (SOI model) with altered death rates
#' @param numberchanges (optional) a vector of the same length as species in the (meta-)community that describes by what amount the abundance of each species is decreased or increased (zero: no change, minus: decrease, else increase, NA: remove a species from the community); species that have reached zero will not be altered
#' @param capacityConstant determines whether the total abundance is allowed to change or not (if capacityConstant is TRUE, the abundance of all non-increasing species will be decreased by an equal amount to keep the total abundance constant)
#' @return a perturbation object
#' @examples
#' # define three perturbations on a community with 10 species
#' per=perturbation(times=c(1,10,15),durations=rep(1,3),numberchanges=c(10,-10,rep(0,8)))
#'
#' @export

perturbation<-function(times=c(), durations=c(), growthchanges=c(), numberchanges=c(), deathrate=NA, capacityConstant=FALSE){
  if(length(times)!=length(durations)){
    stop("Please define for each perturbation its duration.")
  }
  #if(length(growthchanges)==0 && length(numberchanges)==0 && is.na(deathrate)){
  #  stop("Perturbation has no effect on the community. Please provide either the vector of growth rate changes or the vector of abundance changes or specify another death rate.")
  #}
  if(length(deathrate)>1 && (length(growthchanges)>0 || length(numberchanges)>0)){
    if((length(growthchanges)>0  && (length(deathrate) != length(growthchanges))) ||  (length(numberchanges)>0  && (length(deathrate) != length(numberchanges)))){
      stop("If a vector of death rates is provided, please provide as many death rates as growth changes or number changes!")
    }
  }
  perturbation=list(times,durations, growthchanges, numberchanges, deathrate, capacityConstant)
  names(perturbation)=c("times","durations","growthchanges","numberchanges","deathrate","capacityConstant")
  attr(perturbation, "class") <- "perturbation"
  return(perturbation)
}

# apply the perturbation (function called in simulation code)
# perturb: a perturbation object
# t: the current time point in the simulation
# perturbCounter: a counter for the current perturbation
# durationCounter: a counter for the duration of the current perturbation
# perturbationOn: perturbation is active for current time point
# ori.growthrates: original growth rate values
# abundances: current abundance values
applyPerturbation<-function(perturb, t=NA, perturbCounter=1, durationCounter=1, perturbationOn=FALSE, ori.growthrates=c(), abundances=c()){
  growthrates=ori.growthrates
  #print(paste("perturbation status",perturbationOn))
  #print(paste("perturb counter",perturbCounter))
  #print(paste("Length current duration: ",perturb$durations[perturbCounter]))
  #print(paste("perturb time: ",perturb$times[perturbCounter]))
  #print(paste("time: ",t))
  if(perturbCounter<=length(perturb$times)){
    if(perturb$times[perturbCounter]==t){
      perturbationOn=TRUE
      #print("perturbation on")
    }
    if(perturbationOn==TRUE && durationCounter<=perturb$durations[perturbCounter]){
      if(!is.null(perturb$growthchanges)){
        if(durationCounter==1){
          if(length(perturb$growthchanges)>0){
            growthrates=ori.growthrates+perturb$growthchanges
            #print(paste("updated growth rates: ",growthrates))
          }
          indices.rate.neg=which(growthrates<0)
          if(length(indices.rate.neg)>0){
            growthrates[indices.rate.neg]=0 # set growth rates with negative value to zero
          }
        }
      }
      if(!is.null(perturb$numberchanges)){
        if(length(perturb$numberchanges)>0){
          if(perturb$capacityConstant==TRUE){
            posChange.indices=which(perturb$numberchanges>0)
            # only take non-zero indices into account to avoid resurrection of species
            posChange.indices=setdiff(which(abundances>0),posChange.indices)
            totalPosChange=sum(perturb$numberchanges[posChange.indices])
            # to keep total abundance constant, positive change is accompanied by equal
            # negative change in all non-positively changing species
            negChangeNumber=length(perturb$numberchanges)-length(posChange.indices)
            if(negChangeNumber<=0){
              stop("The total abundance cannot be kept constant, since there are no species that can counter-balance the positive change in abundance.")
            }else{
              negChangeAmount=totalPosChange/negChangeNumber
              negChange.indices=setdiff((1:length(perturb$numberchanges)),posChange.indices)
              perturbationValues=perturb$numberchanges
              perturbationValues=perturbationValues[negChange.indices]-negChangeAmount
              prevTotalAbund=sum(abundances)
              abundances=abundances+perturbationValues
              currentTotalAbund=sum(abundances)
              if(prevTotalAbund!=currentTotalAbund){
                warning("The total abundance changed!")
              }
            }
          }else{
            # only change non-zero abundances to avoid species resurrection
            nonzero.indices=which(abundances>0)
            abundances[nonzero.indices]=abundances[nonzero.indices]+perturb$numberchanges[nonzero.indices]
          }
          #print(abundances)
        }
        indices.na=which(is.na(abundances))
        if(length(indices.na)>0){
          abundances[indices.na]=0 # set indicated species to zero
        }
        indices.neg=which(abundances<0)
        if(length(indices.neg)>0){
          abundances[indices.neg]=0 # set species with negative abundance to zero
        }
      }
      durationCounter=durationCounter+1
      #print(paste("Updated duration counter: ",durationCounter))
      #print(paste("perturbation on? ",perturbationOn))
    }else if(perturbationOn==TRUE && durationCounter>perturb$durations[perturbCounter]){
      # duration over, re-apply old values
      growthrates=ori.growthrates
      perturbCounter=perturbCounter+1
      durationCounter=1
      perturbationOn=FALSE
      #print(paste("Updated perturbation counter: ",perturbCounter))
    }
  }
  res=list(growthrates, abundances, perturbCounter, durationCounter, perturbationOn)
  names(res)=c("growthrates","abundances", "perturbCounter","durationCounter", "perturbationOn")
  return(res)
}
