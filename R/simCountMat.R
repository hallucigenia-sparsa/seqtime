#' @title Simulate a count matrix
#'
#' @description Simulate the output of a sequencing experiment, with taxa as rows and samples as columns.
#'
#' @param N the number of taxa
#' @param samples the number of columns in the matrix
#' @param pi taxon proportion vector of length N (if smaller N, it is generated using the mode and k parameters)
#' @param counts either the total number of counts in each sample or a vector specifying the count number in each sample
#' @param distrib "dm" for Dirichlet-Multinomial distribution or "unif" for uniform distribution
#' @param maxcount maximal count number for any taxon (only for uniform distribution)
#' @param mode how the taxon probability vector for the DM distribution is to be generated (1=perfectly even, 2=sampled from uniform distribution, 3=dominant taxon has probability of 0.95 and all others have equal probability, 4=probabilities are sampled from a Poisson distribution with lambda set to 0.5, 5=using bstick function from vegan, 6=using geometric series with parameter k)
#' @param k evenness parameter of mode 6
#' @param theta overdispersion parameter of the DM distribution
#' @param norm normalize matrix column-wise, such that the entries in each column add to one
#' @param shuffle.samples shuffle each sample
#' @return a count matrix or relative abundance matrix
#' @examples tsplot(simCountMat(N=10,samples=50,mode=6,k=0.5),header="Dirichlet-Multinomial")
#' @seealso \code{\link{generateAbundances}}
#' @export

simCountMat<-function(N, pi=c(), samples=100, counts=1000, distrib="dm", maxcount=100, mode=1, k=0.5, theta=0.002, norm=F, shuffle.samples=F){
  if(length(counts)==1 && counts < 1){
    stop("counts should be at least 1.")
  }
  if(N < 1){
    stop("There should be at least one taxon.")
  }
  if(samples < 1){
    stop("There should be at least one sample.")
  }
  if(length(pi) < N){
    pi = generateAbundances(N=N, mode=mode, k=k, probabs=TRUE)
  }
  mat=matrix(nrow=N,ncol=samples)
  # get taxon probabilities for each column
  sample.count = counts
  for(i in 1:samples){
    if(length(counts)==samples){
      sample.count = counts[i]
    }
    mat[,i]=simCounts(pi=pi, theta=theta, counts = sample.count, distrib=distrib, maxcount=maxcount)
    if(shuffle.samples){
      mat[,i]=sample(mat[,i])
    }
  }
  # normalize matrix column-wise
  if(norm){
    mat = normalize(mat)
  }
  mat
}

##############################################################################
# Generate counts with Dirichlet-multinomial (dm) or uniform distribution (unif)
#
# Arguments:
# samples = number of samples to generate
# counts = total number of counts in each sample
# pi = taxon probability vector, its length indicates the number of taxa
# theta = overdispersion parameter for dm
# maxcount = maximal count number for any taxon (only for uniform distribution)
#
# Value:
# a count vector or matrix
#
# Note: taken partly from package HMP, function Dirichlet.multinomial and
# package dirmult, function simPop
#
##############################################################################
simCounts<-function(samples=1,counts=1000,pi=rep(1/10,10),theta=0.002, distrib="dm", maxcount=100){
  if(sum(pi) != 1){
    stop("Taxon probabilities should sum to one!")
  }
  if(theta > 1 || theta < 0){
    stop("Theta should be between 0 and 1.")
  }
  if(distrib == "dm"){
    # from a line in function simPop in package dirmult
    gamma = pi*(1 - theta)/theta
    # from a line in Dirichlet.multinomial from package HMP
    res=rmultinom(n=samples, size=counts, prob=gtools::rdirichlet(samples,gamma))
  }else if(distrib == "unif"){
    if(maxcount < 1){
      stop("The maximal count should be at least 1.")
    }
    res=matrix(runif(length(pi)*samples,0,maxcount),nrow=length(pi),ncol=samples)
    res = round(res)
  }else{
    stop("Choose dm or unif as generating distribution.")
  }
  res
}

