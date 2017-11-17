#' @title Self-organized instable model
#'
#' @description Generate time series with the self-organized instable (SOI) model
#' implementing model B by Sole et al. 2002.
#'
#' @param N number of species
#' @param I number of individuals
#' @param A interaction matrix
#' @param m.vector species-specific immigration probabilities (these also determine initial abundances)
#' @param e.vector species-specific extinction probabilities
#' @param tend number of time points (i.e. the number of generations)
#' @param K TODO
#' @param perturb a perturbation object
#' @return a matrix with species abundances as rows and time points as columns
#' @seealso \code{\link{ricker}} for the Ricker model
#' @references Sole et al. Philos Trans R Soc Lond B Biol. Sci. "Self-organized instability in complex ecosystems" 357:667-671 (2002)
#' @examples
#' \dontrun{
#' N=10
#' A=generateA(N,c=0.1)
#' tsplot(soi(N=N,I=1000,A=A,tend=100), header="SOI")
#' }
#' @export



soi<-function(N, I, A, m.vector=runif(N), e.vector=runif(N), tend, K=5, perturb=NULL){

  results<-generate.parameters(N, I, A, m.vector, e.vector)

  TS<-tend
  omega <-results$A
  speciesNr <-results$N
  I <-results$I
  x0 <- results$abundances
  name_species <-names(x0)
  B <- 0
  for (i in 1:speciesNr){
    B <-B+x0[[i]]
  }
  # species speciesNr+1 corresponds to empty slots
  S <-speciesNr+1
  x0 <-c(x0, I-B)
  names(x0) <-c(name_species, paste("Spec_",S,"",sep=""))


  # now we prepare the call to run.SOI with the current parameters,
  # we compute a list of jumps size and the parameters for the propensity computation
  nu <- matrix(0*(1:(2*S*speciesNr)), nrow=S)
  help_calc=list()
  coeff_trans <-c()
  name_coeff <- c()
  helpful<-c()
  names(coeff_trans)<-name_coeff
  # we build the transisiton matrix
  # immigration & positive interaction jumps
  for (i in 1:speciesNr){
    nu[i,i]<-1
    nu[S,i]<- -1
    weaker<-which((omega[i,]>omega[,i]) & omega[,i]>=0)
    ##fix:
    if (length(weaker)==0){
      c_prop<-integer(0)
    }
    else {
      c_prop<-omega[i,weaker]+omega[weaker,i]
    }
    #c_prop<-omega[i,weaker]+omega[weaker,i]
    helpful<-rbind(weaker,c_prop)
    help_calc[[i]]<-helpful
  }
  # extinction jumps
  for (i in 1:speciesNr){
    nu[i,(i+speciesNr)]<- -1
    nu[S,(i+speciesNr)]<- 1
    help_calc[[speciesNr+i]]=integer(0)
  }
  # interaction jumps

  for (i in 1:speciesNr){
    for (j in 1:speciesNr){
      if(omega[j,i]<0 & omega[j,i]<omega[i,j]){
        competition<-1
        pji=omega[i,j]-omega[j,i]
        helpful<-c(i,j,competition, pji)
        help_calc<-c(help_calc, list(helpful))
        jump <- 0*(1:S)
        jump[i]<- 1
        jump[j]<- -1
        nu<-cbind(nu,jump)}
    }}
  # now we build parms (needed to run rateFunc, which computes the propensity):
  #I, S, mu, e, help_calc
  parms=list()
  parms[[1]]=I
  parms[[2]]=speciesNr
  parms[[3]]=results$immigration_prob
  parms[[4]]=results$extinction_prob
  parms[[5]]=help_calc

  #call to the simulation
  #the simulation method is the K-leap method, no test on K.
  abundances_over_time<-run.SOI(TS,x0,nu,parms,K,perturb)
  #output provides the simulated abundances (colums) at each time (rows), excluding the
  #initial time
  return(t(abundances_over_time))

}



###############
generate.parameters<-function(N, I, A, m.vector, e.vector){
  # species_list
  species_list <- list()
  for (i in 1:N) {
    species_list[[paste("Spec_",i,sep="")]]<-0
  }

  n_spec<-c()
  for (i in 1:N) {
    n_spec<-c(n_spec, paste("Spec_", i, sep = ""))
  }
  names(m.vector)<-n_spec
  names(e.vector)<-n_spec

  # start with partially filled sited list
  # sites are filled with randomly picked species
  # according to their immigration probabilities
  sites<-c()
  n<-c()
  for (i in 1:I){
    species<-sample(1:N,1)
    if (runif(1,min=0,max=1)<m.vector[species]){
      sites[i]<-species
    } else sites[i]<-0
    n<-c(n,paste("site_",i,sep=""))
  }
  names(sites)<-n

  # generate abundance list (species vs. abundance)
  abundances<-c()
  n_spec<-c()
  for (i in 1:N) {
    abundances[i]<-sum(sites == i)
    n_spec<-c(n_spec, paste("Spec_", i, sep = ""))
  }
  names(abundances)<-n_spec

  # generate/update species_list and living_species
  living_species<-vector()
  for (i in 1:N){
    occurrance<-which(sites == i)
    species_list[[i]]<-occurrance
    if (abundances[[i]]!=0){
      living_species<-append(living_species,i)
    }
  }

  parameters <- list(A=A,species_list=species_list,
                     abundances=abundances,immigration_prob=m.vector,
                     extinction_prob=e.vector,sites=sites,
                     living_species=living_species,I=I,N=N)
}



##########################
rateFunc <- function(x, params, t) {
  # params[[1]]: I number of individuals
  # params[[2]]: S number of actual species
  # params[[3]]: mu immigration rates
  # params[[4]]: e extinction rates
  # params[[5]]: help_calc stuff needed to compute interactions
  propensity <- c()
  S=params[[2]]
  I=params[[1]]
  # we fill the propensity vector, same order as the reaction jump matrix
  # first immigration and mutualistic interactions
  propensity<-c(propensity, x[S+1]*params[[3]]/S )
  for (i in 1:S){
    M=params[[5]][[i]]
    if (ncol(M)>0) {
      propensity[i]=propensity[i]+sum(x[M[1,]]*M[2,])/I*x[i]/I*x[S+1]
    }
  }

  # then extinction
  propensity<-c(propensity, x[1:S]*params[[4]])
  # and finally competition
  for (l in (2*S+1):length(params[[5]])){
    i<-params[[5]][[l]][1]
    j<-params[[5]][[l]][2]
    pji<-params[[5]][[l]][4]
    propensity<-c(propensity,x[i]*x[j]*pji/I)
  }
  names(propensity)<-NULL
  return(propensity)
}


#########################
run.SOI<-function(time_steps,x0,nu,parms,K,perturb){
  t<-0
  ######
  abundances_over_time <- matrix(nrow=time_steps, ncol=parms[[2]])
  ######
  x<-x0
  index<-1
  #############
  while(index<=time_steps){
    propensity<-rateFunc(x,parms,t)
    a0<-sum(propensity)
    theta=propensity/a0
    Km=rmultinom(1,K,theta)
    tau<-rgamma(1,shape=K,rate=a0)
    y<-x+nu%*%Km
    z<-which(y<0)
    n<- -sum(y[z])
    if(n>0){
      y[z]<-0*z
      for (j in 1:n){
        m<-which(y>=1)
        u<-sample(m,1)
        y[u]<-y[u]-1
      }}
    t<-t+tau
    if (t>=index){
      z<-t(y)
      abundances_over_time[index,]<-z[1:parms[[2]]]
      index<-index+1
    }
    x<-y
  }
  return(abundances_over_time)
}

