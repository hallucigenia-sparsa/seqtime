#' Adjust connectance
#'
#' The connectance is adjusted by randomly adding or removing interactions.
#' @param A interaction matrix
#' @param c target connectance
#' @return adjusted interaction matrix
#' @examples
#' A=generateA(N=10)
#' A.adj=adjustc(A,c=0.5)
#' @export

adjustc<-function(A,c=0.02){
  c_obs = getConnectance(A)
  print(paste("Initial edge number", length(A[A!=0])))
  print(paste("Initial connectance", c_obs))
  if(c_obs < c){
    edgeNumAdded = 0
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
        A[xpos,ypos]=1
        c_obs=getConnectance(A=A,N=N)
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
  return(A)
}
