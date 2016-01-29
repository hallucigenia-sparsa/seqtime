#' Modify the interaction matrix
#'
#' Mode negpercent: The number of interactions equals the number of non-zero entries in the interaction matrix.
#' If symmetric is true, both matrix triangles will be set to a negative value (corresponding
#' to competition), else only one triangle will be set to a negative value (corresponding
#' to parasitism, predation or amensalism). If symmetric is true, the number of interactions
#' is counted only in one triangle of the interaction matrix.
#'
#' @param A the interaction matrix
#' @param mode modification mode, adjustc (adjust connectance to reach target connectance), negpercent (set the specified percentage of negative edges), tweakstable (randomly add negative edges until matrix is stable or fully connected)
#' @param strength interaction strength, binary (0/1) or uniform (sampled from uniform distribution from minstrength to 1)
#' @param minstrength minimum interaction strength for uniform mode (maximum is 1)
#' @param c the target connectance (only for mode adjustc)
#' @param perc negative edge percentage (only for mode negpercent)
#' @param symmetric only introduce symmetric negative interactions (only for mode negpercent)
#' @return the modified interaction matrix
#'
modifyA<-function(A, mode="adjustc", strength="binary", minstrength=0.1, c=0.2, perc=50, symmetric=FALSE){
  edgeNumAdded = 0
  print(paste("Initial edge number", length(A[A!=0])))
  c_obs = getConnectance(A)
  print(paste("Initial connectance", c_obs))
  if(mode == "adjustc"){
    if(c_obs < c){
      while(c_obs < c){
        # randomly select source node of edge
        xpos=sample(c(1:ncol(A)))[1]
        # randomly select target node of edge
        ypos=sample(c(1:ncol(A)))[1]
        # avoid diagonal
        if(xpos != ypos){
          # count as added if there was no edge yet
          if(A[xpos,ypos]==0){
            edgeNumAdded = edgeNumAdded+1
          }
          # add edge
          A[xpos,ypos]=getStrength(strength=strength,pos=TRUE, minstrength=minstrength)
          c_obs=getConnectance(A=A)
        }
      }
      print(paste("Number of edges added", edgeNumAdded))
    }else if(c_obs > c){
      edgeNumRemoved = 0
      while(c_obs > c){
        xpos=sample(c(1:ncol(A)))[1]
        ypos=sample(c(1:ncol(A)))[1]
        # avoid diagonal
        if(xpos != ypos){
          # count as removed if there was an edge before
          if(A[xpos,ypos]!=0){
            edgeNumRemoved = edgeNumRemoved+1
          }
          # remove edge
          A[xpos,ypos]=0
          c_obs = getConnectance(A)
        }
      }
      print(paste("Number of edges removed", edgeNumRemoved))
    }
    print(paste("Final connectance", c_obs))
  }else if(mode == "negpercent"){
    # arc number: number of non-zero entries in the interaction matrix
    num.edge=length(A[A!=0])
    num.perc=(num.edge/100)*perc
    # symmetric interactions: we will count negative edges, not arcs
    if(symmetric == TRUE){
      num.perc=round(num.perc/2)
    }
    print(paste("Converting",num.perc,"edges into negative edges",sep=" "))
    indices=which(A==1,arr.ind=T)
    # randomly select indices
    rand=sample(c(1:nrow(indices)))
    xyposAlreadySeen = c()
    counter = 0
    # loop over number of negative edges to introduce
    for(i in 1:num.perc){
      xpos=indices[rand[i],1]
      ypos=indices[rand[i],2]
      xypos=paste(xpos,"_",ypos, sep="")
      yxpos=paste(ypos,"_",xpos,sep="")
      # if we find an index pair that was already used, we have to look for another index pair,
      # since using the same index pair means to use the same arc or the same arc in reverse direction
      if(symmetric == TRUE && is.element(xypos,xyposAlreadySeen) == TRUE){
        xpos = indices[rand[nrow(indices)-counter],1]
        ypos = indices[rand[nrow(indices)-counter],2]
        counter = counter + 1
        if((num.perc + counter) > nrow(indices)){
          stop("More negative edges requested than can be set!")
        }
      }
      xyposAlreadySeen = c(xypos, yxpos, xyposAlreadySeen)
      # print for tests
      # print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
      negval=getStrength(strength=strength,pos=FALSE,minstrength=minstrength)
      A[xpos,ypos]=negval
      if(symmetric == TRUE){
        A[ypos,xpos]=negval
      }
      #print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
      #print(paste("reverse value:",A[ypos,xpos],sep=" "))
    }
  }else if(mode == "tweakstable"){
    # add negative edges randomly until matrix is stable or fully connected
    while(testStability(A, method="ricker")==FALSE && isFullyconnected(A)==FALSE){
      # randomly select source node of edge
      xpos=sample(c(1:ncol(A)))[1]
      # randomly select target node of edge
      ypos=sample(c(1:ncol(A)))[1]
      # avoid diagonal
      if(xpos != ypos){
        # count as added if there was no edge yet
        if(A[xpos,ypos]==0){
          edgeNumAdded = edgeNumAdded+1
        }
        # add negative arc with random interaction strength
        A[xpos,ypos]=getStrength(strength=strength,pos=FALSE,minstrength=minstrength)
        #print(paste("Number of edges added", edgeNumAdded))
      }
    } # end tweak A
    print(paste("Number of edges added", edgeNumAdded))
  }
  c=getConnectance(A)
  print(paste("Final connectance:",c))
  return(A)
}

################## helper functions ################

# Get the interaction strength.
getStrength<-function(strength="binary", minstrength=0.1, pos=TRUE){
  value = NA
  if(strength=="binary"){
    value = 1
  }else if(strength == "uniform"){
    value = runif(1,min=minstrength,max=1)
  }
  if(!pos){
    value = -1*value
  }
  return(value)
}

# Check whether the interaction matrix is fully connected (entirely filled with 1)
isFullyconnected<-function(A){
  if(length(A[A!=0])==nrow(A)*nrow(A)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
