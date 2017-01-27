#' Modify the interaction matrix
#'
#' Mode negpercent: The number of interactions equals the number of non-zero entries in the interaction matrix.
#' If symmetric is true, both matrix triangles will be set to a negative value (corresponding
#' to competition), else only one triangle will be set to a negative value (corresponding
#' to parasitism, predation or amensalism). If symmetric is true, the number of interactions
#' is counted only in one triangle of the interaction matrix.
#'
#' @param A the interaction matrix
#' @param mode modification mode, values: adjustc (adjust connectance to reach target connectance), schur (remove positive eigenvalues using the Schur decomposition), negpercent (set the specified percentage of negative edges), tweak (multiply a randomly chosen positive interaction strength with -1), enforceneg (multiply negative interaction strengths with given factor, but keep diagonal as is), removeorphans (remove all species that do not interact with other species), mergeposlinks, mergeneglinks, mergelinks (merge positive/negative/all links of all taxa with the same name)
#' @param strength interaction strength, binary (0/1) or uniform (sampled from uniform distribution from minstrength to 1)
#' @param factor multiplication factor for enforceneg mode
#' @param minstrength minimum interaction strength for uniform mode (maximum is 1)
#' @param c the target connectance (only for mode adjustc)
#' @param perc negative edge percentage (only for mode negpercent)
#' @param symmetric only introduce symmetric negative interactions (only for mode negpercent)
#' @return the modified interaction matrix
#' @export

modifyA<-function(A, mode="adjustc", strength="binary", factor=2, minstrength=0.1, c=0.2, perc=50, symmetric=FALSE){
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
  }else if(mode=="mergeposlinks" || mode=="mergeneglinks" || mode=="mergelinks"){
    entries=unique(rownames(A))
    mergedlinks=matrix(0,nrow=length(entries),ncol=length(entries))
    rownames(mergedlinks)=entries
    colnames(mergedlinks)=entries
    for(i in 1 : nrow(A)){
      for(j in 1 : ncol(A)){
        if(!is.na(A[i,j]) && A[i,j]!=0){
          merge=FALSE
          xIndex=which(entries==rownames(A)[i])
          yIndex=which(entries==rownames(A)[j])
          #print(paste("x index:",xIndex))
          #print(paste("y index:",yIndex))
          if(mode=="mergeposlinks" && A[i,j]>0){
            merge=TRUE
          }else if(mode=="mergeneglinks" && A[i,j]<0){
            merge=TRUE
          }else if(mode=="mergelinks"){
            merge=TRUE
          }
          if(merge==TRUE){
            mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]+1
          }
        } # interaction is not zero
      } # column loop
    } # row loop
    A=mergedlinks
  }else if(mode=="removeorphans"){
    # since A can be asymmetric, only those species can be removed for which rows and columns are simultaneously zero (except for diagonal)
    toKeep=c()
    diagvals=diag(A)
    diag(A)=0
    for(i in 1:nrow(A)){
      rowsum=sum(abs(A[i,]))
      colsum=sum(abs(A[,i]))
      if(rowsum != 0 || colsum!=0){
        toKeep=append(toKeep,i)
      }
    }
    A=A[toKeep,toKeep]
    diag(A)=diagvals[toKeep]
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
  }else if(mode == "tweak"){
    # check that positive arcs are present
    if(length(A[A>0]) > 0){
      # row and column indices of positive arcs
      indices.pos = which(A>0,arr.ind=TRUE)
      # randomly select a positive arc
      x=sample(1:nrow(indices.pos),1)
      # convert positive arc into negative one, keeping the same interaction strength
      A[indices.pos[x,1],indices.pos[x,2]]=A[indices.pos[x,1],indices.pos[x,2]]*(-1)
    }else{
      warning("Cannot tweak. No positive arc in the given matrix.")
    }
  }else if(mode == "enforceneg"){
    diag=diag(A)
    indices.neg = which(A<0,arr.ind=TRUE)
    # multiply negative entries by given factor
    A[indices.neg]=A[indices.neg]*factor
    # keep original diagonal
    diag(A)=diag
  }else if(mode == "schur"){
    # remove positive real parts of eigenvalues if any (using schur decomposition)
    sd<-dim(A)

    if(max(Re(eigen(A)$values))){
      # division by max.A helps removing positive eigenvalues
      max=max(A)
      A=A/max

      diagt<-diag(sd[2])+0i

      # Computes the generalized eigenvalues and Schur form of a pair of matrices.
      # R=imaginary part identical to 0 with a tolerance of 100*machine_precision as determined by Lapack
      schur.A<-geigen::gqz(A,diagt,"R")
      # generalized inverse of a matrix
      T<-schur.A$S%*%MASS::ginv(schur.A$T)
      rediag<-Re(diag(T))
      imdiag<-Im(diag(T))

      indicesP=rediag>0
      listind=1:sd[2]

      for(k in listind[indicesP]){
        T[k,k]<- complex(real=-Re(T[k,k]),imaginary=Im(T[k,k]))
      }

      A <- schur.A$Q %*% T %*% MASS::ginv(schur.A$Q)
      A<-Re(A)
      A=A*max
    }

  }else{
    stop(paste("Mode",mode,"not known."))
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
