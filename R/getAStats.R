#' @title Analyse an interaction matrix
#'
#' @description Count interaction types or compute network properties of interaction matrix.
#' The mean interaction strength is computed according to Coyte and colleagues
#' by fitting a half-normal distribution to the realized interaction strengths.
#' Graph properties (modularity, average clustering coefficient, average path length)
#' are computed using igraph functions.
#'
#' @param A the interaction matrix
#' @param statsType interactions, degree or network
#' @param plot.degree plot the degree distribution (for statsType network)
#' @param collapse.degree sum degrees for taxa with the same name (for statsType degree)
#' @return for degree a matrix with positive, negative and total degree (including self-loops, excluding missing values), else a list; for statsType interactions: meanstrength = average interaction strength, varstrength = variance of interaction strength, nbinteractions = total interaction number (excluding diagonal), nbmut = number of mutualisms, nbcomp = number of competitions, nbcom = number of commensalisms, nbam = number of amensalisms, nbexp = number of exploitations, for statsType network: nodenum (node number), arcnum (arc number, including diagonal), mod (fast greedy modularity), cc (average clustering coefficient), avgpathlength (average shortest path length)
#' @references Coyte et al., Science: "The ecology of the microbiome: Networks, competition, and stability" 350 (6261), 663-666 (2015).
#' @export

getAStats<-function(A, statsType = "interactions", plot.degree=FALSE, collapse.degree=FALSE){

  if(statsType == "interactions"){
    # numbers of interaction types, excluding diagonal
    P=0
    Pe=0 # +/- exploitative (predator/prey, parasite/host)
    Pc=0 # -/- competition
    Pm=0 # +/+ mutualism
    Pp=0 # +/0 commensalism
    Pa=0 # -/0 amensalism
    # non-paired-zero, non-diagonal interaction values
    a = c()

    # collect values
    for(i in 1:ncol(A)){
      # skip upper triangle
      for(j in 1:i){
        # skip diagonal
        if(i != j){
          #print(paste("i=",i))
          #print(paste("j=",j))
          # collect non-zero interaction strengths
          if((A[i,j] != 0) || (A[j,i] != 0)){
            a = c(a,A[i,j])
          }
          # exploitative
          if((A[i,j] > 0 && A[j,i] < 0) || (A[i,j] < 0 && A[j,i] > 0)){
            Pe = Pe + 1
          }
          # competitive (symmetric)
          else if(A[i,j] < 0 && A[j,i] < 0){
            Pc = Pc + 1
          }
          # mutualistic (symmetric)
          else if(A[i,j] > 0 && A[j,i] > 0){
            Pm = Pm + 1
          }
          # commensalistic
          else if((A[i,j] > 0 && A[j,i] == 0) || (A[i,j] == 0 && A[j,i] > 0)){
            Pp = Pp + 1
          }
          # amensalistic
          else if((A[i,j] < 0 && A[j,i] == 0) || (A[i,j] == 0 && A[j,i] < 0)){
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

    P=Pm+Pc+Pe+Pp+Pa

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
  else if(statsType=="degree"){
    degrees=matrix(nrow=nrow(A),ncol=3)
    A[is.na(A)]=0
    # includes loops
    for(i in 1:nrow(A)){
      diagVal=A[i,i]
      # outgoing plus incoming arcs are counted
      degrees[i,1]=length(which(A[i,]>0))+length(which(A[,i]>0))
      if(diagVal>0 && degrees[i,1]!=0){
        # subtract diagonal, else diagonal is counted twice (for incoming and outgoing arcs)
        degrees[i,1]=degrees[i,1]-1
      }
      # outgoing plus incoming arcs are counted
      degrees[i,2]=length(which(A[i,]<0))+length(which(A[,i]<0))
      if(diagVal < 0 && degrees[i,2]!=0){
        degrees[i,2]=degrees[i,2]-1
      }
      degrees[i,3]=degrees[i,1]+degrees[i,2]
    }
    rownames(degrees)=rownames(A)
    posindex=which(degrees[,1]==max(degrees[,1]))
    print(paste(rownames(degrees)[posindex]," has a maximal positive degree of ",max(degrees[,1]),sep=""))
    negindex=which(degrees[,2]==max(degrees[,2]))
    print(paste(rownames(degrees)[negindex]," has a maximal negative degree of ",max(degrees[,2]),sep=""))
    index=which(degrees[,3]==max(degrees[,3]))
    print(paste(rownames(degrees)[index]," has a maximal total degree of ",max(degrees[,3]),sep=""))

    if(collapse.degree==TRUE){
      entries=unique(rownames(degrees))
      collapsed=matrix(nrow=length(entries),ncol=ncol(degrees))
      entryCounter=1
      for(entry in entries){
        print(paste("Processing entry",entry))
        indices=which(rownames(degrees)==entry)
        for(colIndex in 1:ncol(degrees)){
          collapsed[entryCounter,colIndex]=sum(degrees[indices,colIndex])
        }
        entryCounter=entryCounter+1
      }
      rownames(collapsed)=entries
      degrees=collapsed
    }
    colnames(degrees)=c("pos","neg","total")
    res=degrees
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
      plot(degree.distribution(g,cumulative=T),log="xy")
    }
    cc=transitivity(g, type="average")
    avglength=average.path.length(g)
    print(paste("Average clustering coefficient of final network:",cc))
    print(paste("Average path length of final network:",avglength))

    res=list(nodenum,arcnum,maxFc,cc,avglength)
    names(res)=c("nodenum","arcnum","mod","cc","avgpathlength")

  }
  return(res)
}
