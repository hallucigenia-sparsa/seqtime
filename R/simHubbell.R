#' Simulate species abundances with the basic Hubbell model.
#'
#' @param N species number
#' @param I number of individuals
#' @param y initial species probabilities, should sum to one
#' @param m.vector species proportions in the metacommunity, should sum to one
#' @param m immigration rate (probability that a dead individuum is replaced by an individuum from the metacommunity)
#' @param d number of deaths at each time step
#' @param tend number of time points (i.e. the number of generations)
#' @return a matrix with species abundances as rows and time points as columns
#' @example tsplot(simHubbell(N=100,I=3000,d=300))
#' @export

simHubbell<-function(N=50, I=500, y=rep(1/N,N), m.vector=rep(1/N,N), m=0.02, d=10, tend=100){
  if(sum(m.vector) != 1){
    stop("Species proportions in the metacommunity need to sum to one.")
  }
  if(sum(y) != 1){
    stop("Initial species probabilities need to sum to one.")
  }

  gridsize=round(sqrt(I))
  # fill grid with ones
  grid=matrix(1,nrow=gridsize, ncol=gridsize)
  S=c(1:N)
  # initialize the grid such that each species is selected according to its initial probability
  for(i in 1:gridsize){
    for(j in 1:gridsize){
      grid[i,j]=sample(S,1,prob=y)
    }
  }

  counts=countSpecies(grid,N)

  tseries=matrix(NA,nrow=N,ncol=tend)
  tseries[,1]=counts

  for(t in 2:tend){
    z1=ceiling(gridsize*runif(d)) # x positions in the grid, ceil is like round, but guarantees no integer smaller than 1 will occur
    z2=ceiling(gridsize*runif(d)) # y positions in the grid
    z=runif(d) # generate d uniformly distributed values
    immigrants=as.numeric(z>=m)
    locals=as.numeric(z<m)
    # immigration probabs
    immiprobabs=immigrants*randp(d,m.vector,S)
    # local replacement probabs
    localprobabs=locals*randp(d,counts,S)
    r=immiprobabs+localprobabs
    # replace dead individuals with local or immigrated ones
    for(i in 1:d){
      grid[z1[i],z2[i]]=r[i]
    }

    # update counts
    counts=countSpecies(grid,N)

    tseries[,t]=counts

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

# if indplot is true, plot the number of individuals over time,
# else plot the number of species over time
# example: diagnostic(simHubbell(N=100,tend=500),indplot=TRUE)
# community size: species number
# for simHubbell:
# number of individuals fluctuates, number of species always goes down for the first 1000 time points
# for simUNTB:
# number of individuals stays constant, number of species tends to increase for the first 1000 time points, but not always
diagnostic<-function(x, indplot=TRUE){
  if(indplot == TRUE){
    colsums=apply(x,2,sum)
    plot(1:ncol(x),colsums, ylab="Number of individuals", xlab="Time points") # no trend?
  }else{
    # species number bias - species number decreases with time
    specnum=c()
    for(j in 1:ncol(x)){
      v=as.numeric(x[,j])
      specnum=c(specnum, length(v[v>0]))
    }
    plot(1:ncol(x),specnum, ylab="Number of species", xlab="Time points")
  }
}
