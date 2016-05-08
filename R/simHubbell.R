#' Simulate species abundances with the basic Hubbell model.
#'
#' For a recent review on the Hubbell model, see Rosindell, Hubbell and Etienne,
#' The Unified Neutral Theory of Biodiversity and Biogeography at Age Ten,
#' Trends in Ecology and Evolution, vol. 26 (7), 340-348, 2011.
#'
#' @param N species number in the local community
#' @param M species number in the metacommunity
#' @param I number of individuals
#' @param y initial species probabilities, should sum to one
#' @param m.vector species proportions in the metacommunity, should sum to one
#' @param m immigration rate (probability that a dead individuum is replaced by an individuum from the metacommunity)
#' @param d number of deaths at each time step
#' @param tskip number of initial time points to be skipped when returning the result (to avoid the transient)
#' @param tend number of time points (i.e. the number of generations)
#' @return a matrix with species abundances as rows and time points as columns
#' @example
#' N=50
#' M=500
#' metapop=generateAbundances(N=M, mode=5, probabs=TRUE)
#' tsplot(simHubbell(N=N, M=M,I=3000,d=N, m.vector=metapop, tskip=500, tend=1000))
#' @export

simHubbell<-function(N=50, M=500, I=500, y=rep(1/N,N), m.vector=rep(1/M,M), m=0.02, d=10, tskip=0, tend=100){
  if(sum(m.vector) != 1){
    stop("Species proportions in the metacommunity need to sum to one.")
  }
  if(sum(y) != 1){
    stop("Initial species probabilities need to sum to one.")
  }
  if(tskip >= tend){
    stop("The period to be skipped is as large as the entire simulation period!")
  }

  gridsize=round(sqrt(I))
  # fill grid with ones
  grid=matrix(1,nrow=gridsize, ncol=gridsize)
  S=c(1:N)
  Smeta=c(1:M)
  # initialize the grid such that each species is selected according to its initial probability
  for(i in 1:gridsize){
    for(j in 1:gridsize){
      grid[i,j]=sample(S,1,prob=y)
    }
  }

  # count up to M, since new species may appear later from the metacommunity
  counts=countSpecies(grid,M)

  tseries=matrix(NA,nrow=M,ncol=(tend-tskip))
  tseriesCounter=1

  if(tskip==0){
    tseries[,tseriesCounter]=counts
    tseriesCounter=tseriesCounter+1
  }

  for(t in 2:tend){
    z1=ceiling(gridsize*runif(d)) # x positions in the grid, ceil is like round, but guarantees no integer smaller than 1 will occur
    z2=ceiling(gridsize*runif(d)) # y positions in the grid
    z=runif(d) # generate d uniformly distributed values
    immigrants=as.numeric(z>=m)
    locals=as.numeric(z<m)
    # immigration probabs
    immiprobabs=immigrants*randp(d,m.vector,Smeta)
    # local replacement probabs
    localprobabs=locals*randp(d,(counts/sum(counts)),Smeta)
    r=immiprobabs+localprobabs
    # replace dead individuals with local or immigrated ones
    for(i in 1:d){
      grid[z1[i],z2[i]]=r[i]
    }

    # update counts
    counts=countSpecies(grid,M)

    if(t > tskip){
      tseries[,tseriesCounter]=counts
      tseriesCounter=tseriesCounter+1
    }

  } # simulation until end
  return(tseries)
}

# count the number of occurrences of each species
# species are indexed from 1 to N
countSpecies<-function(grid,N){
  counts=c()
  for(i in 1:N){
    indices=which(grid==i)
    count=0
    if(length(indices) > 0){
      count=length(indices)
    }
    counts=c(counts,count)
  }
  return(counts)
}

# Generate d integer values (d: number of dead individuals)
# where integers from vector v are selected with
# probability m.vector.
randp<-function(d,m.vector,v){
  if(length(m.vector) != length(v)){
    stop("randp: Both vectors must have the same length!")
  }
  r=c()
  for(i in 1:d){
    r=c(r,sample(v,1,prob=m.vector))
  }
  return(r)
}

