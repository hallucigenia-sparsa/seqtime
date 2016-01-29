#' Analyse an interaction matrix
#'
#' Count interaction types or compute network properties of interaction matrix.
#' The mean interaction strength is computed according to Coyte et al.
#' by fitting a half-normal distribution to the realized interaction strengths.
#' Graph properties (modularity, average clustering coefficient, average path length)
#' are computed using igraph functions.
#'
#' @param A the interaction matrix
#' @param statsType interactions or network
#' @param plot.degree plot the degree distribution (for statsType network)
#' @return a list, for statsType interactions: meanstrength = average interaction strength, varstrength = variance of interaction strength, nbinteractions = total interaction number,nbmut = number of mutualisms, nbcomp = number of competitions, nbcom = number of commensalisms, nbam = number of amensalisms, nbexp = number of exploitations, for statsType network: nodenum (node number), arcnum (arc number), mod (fast greedy modularity), cc (average clustering coefficient), avgpathlength (average shortest path length)
#' @references Coyte et al., Science 2015: "The ecology of the microbiome: Networks, competition, and stability"
#' @export

getAStats<-function(A, statsType = "interactions", plot.degree=FALSE){

  if(statsType == "interactions"){
    # proportions of interaction types, excluding diagonal
    Pe=0 # +/- exploitative (predator/prey, parasite/host)
    Pc=0 # -/- competition
    Pm=0 # +/+ mutualism
    Pp=0 # +/0 commensalism
    Pa=0 # -/0 amensalism
    # non-paired-zero, non-diagonal interaction values
    a = c()

    # collect values
    for(i in 1:nrow(A)){
      # skip upper triangle
      for(j in 1:i){
        # skip diagonal
        if(i != j){
          # collect non-zero interaction strengths
          if((A[i,j] != 0) || (A[j,i] != 0)){
            a = c(a,A[i,j])
          }
          # exploitative
          if((A[i,j] > 0 && A[j,i] < 0) || (A[i,j] < 0 && A[j,i] > 0)){
            Pe = Pe + 1
          }
          # competitive (symmetric)
          if(A[i,j] < 0 && A[j,i] < 0){
            Pc = Pc + 1
          }
          # mutualistic (symmetric)
          if(A[i,j] > 0 && A[j,i] > 0){
            Pm = Pm + 1
          }
          # commensalistic
          if((A[i,j] > 0 && A[j,i] == 0) || (A[i,j] == 0 && A[j,i] > 0)){
            Pp = Pp + 1
          }
          # amensalistic
          if((A[i,j] < 0 && A[j,i] == 0) || (A[i,j] == 0 && A[j,i] < 0)){
            Pa = Pa + 1
          }
        } # i != j
      } # end second loop
    } # end first loop

    # compute average interaction strength
    a=abs(a)
    var=var(a)
    # mean strength of realized interactions
    EX=sqrt((2*var)/pi)
    EXsquare=(2*var)/pi

    P=Pa+Pc+Pe+Pm+Pp

    if(length(a) != P){
      stop("The number of interactions should equal the sum of the interaction type numbers!")
    }

    # report statistics
    print(paste("Total number of inter-species interactions:",P))
    print(paste("Number of mutualistic inter-species interactions:",Pm))
    print(paste("Number of competitive inter-species interactions:",Pc))
    print(paste("Number of exploitative inter-species interactions:",Pe))
    print(paste("Number of commensalistic inter-species interactions:",Pp))
    print(paste("Number of amensalistic inter-species interactions:",Pa))
    print(paste("Average inter-species interaction strength:",EX))
    print(paste("Variance of inter-species interaction strengths:",var))

    res=list(EX,var,P,Pm,Pc,Pe,Pp,Pa)
    names(res)=c("meanstrength","varstrength","nbinteractions","nbmut","nbcomp","nbexp","nbcom","nbam")

  }
  else if(statsType == "network"){
    # has to be 0 or 1
    A[A!=0]=1
    g=graph.adjacency(A, mode="directed")
    nodenum=vcount(g)
    arcnum=ecount(g)
    fc <- fastgreedy.community(as.undirected(g),weights=NULL)
    maxFc=max(fc$modularity)
    print(paste("Fastgreedy modularity of final network:",maxFc))
    if(plot.degree==TRUE){
      plot(degree.distribution(g,cum=T),log="xy")
    }
    cc=transitivity(g, type="average")
    avglength=average.path.length(g)
    print(paste("Average clustering coefficient of final network:",cc))
    print(paste("Average path length of final network",avglength))

    res=list(nodenum,edgenum,maxFC,cc,avglength)
    names(res)=c("nodenum","arcnum","mod","cc","avgpathlength")

  }
  return(res)
}
