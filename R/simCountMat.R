#' Simulate a count matrix
#'
#' Simulate the output of a sequencing experiment, with taxa as rows and samples as columns.
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
#' @examples tsplot(simCountMat(taxa=10,samples=50,mode=6,k=0.5),type="l", header="Dirichlet-Multinomial")
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
    pi = getProbabs(taxa=N, mode=mode, k=k)
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

######################## Helper functions ###################################

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

##############################################################################
# Generate a taxon probability vector
#
# Arguments:
# taxa = the length of the taxon probability vector
# mode = 1-6
#        1 = perfect evenness (pi = rep(1/N,N))
#        2 = probabilities sampled from uniform distribution and normalized to one
#        3 = dominant species has a probability of 0.95 and all other species are equally distributed to add probabilities up to 1
#        4 = probabilities are sampled from a poisson distribution with lambda set to 0.5
#        5 = probabilities are generated using broken-stick function bstick() from vegan
#        6 = probabilities are generated with the geometric series of given evenness k
# k    = evenness parameter for mode 6
#
# Value:
# A taxon probability vector
#
##############################################################################
getProbabs<-function(taxa = 10, mode = 1, k = 0.5){
  if(taxa < 1){
    stop("taxa should be at least 1.")
  }
  if(!is.wholenumber(mode)){
    stop("The mode should be an integer between 1 and 6")
  }
  if(mode == 1){
    pi=rep(1/taxa,taxa)
  }else if(mode == 2){
    pi = runif(taxa)
    pi =pi/sum(pi)
  }else if(mode == 3){
    dominant=0.95
    commoner=0.05/(taxa-1)
    pi=c(dominant,rep(commoner,(taxa-1)))
  }else if(mode == 4){
    pi=rpois(taxa,lambda=0.5)
    pi=pi/sum(pi)
  }else if(mode == 5){
    pi = vegan::bstick(taxa,1)
  }else if(mode == 6){
    pi=geom.series(l=taxa, counts=1, k = k)
  }else if(mode < 1 || mode > 6){
    stop("Choose a mode between 1 and 6.")
  }
  pi
}

##############################################################################
# Check whether a number is an integer.
#
# Argument:
# x = number
# tol = tolerance
#
# Value:
# boolean true if number is an integer within the tolerance, false otherwise
#
# Note: function taken from http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
#
##############################################################################
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

##############################################################################
# Generate geometric series (niche apportionment model).
#
# Arguments:
# l = length of geometric series
# counts = sum of geometric series
# k = evenness of geometric series, between 0 and 1 (the smaller, the more even)
#
# Value:
# The geometric series
#
# Note: The geometric series is for instance described in
# "Comments about some species abundance patterns: classic, neutral, and niche partitioning models"
# by Ferreira and Petrere-Jr., Braz. J. Biol. 2008
#
##############################################################################
geom.series<-function(l=10, counts=1000, k=0.5){
  if(l < 1){
    stop("l should be at least 1.")
  }
  if(k > 1 || k < 0){
    stop("k should be a number between 0 and 1.")
  }
  pi=c()
  C = 1/(1 - ((1-k)^l))
  for(i in 1:l){
    pi=c(pi,counts*C*k*(1-k)^(i-1))
  }
  pi=pi/sum(pi)
  pi
}
