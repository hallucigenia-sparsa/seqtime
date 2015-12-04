#' Set given percentage of interactions to negative values.
#'
#' The number of interactions equals the number of non-zero entries in the interaction matrix.
#' If symmetric is true, both matrix triangles will be set to a negative value (corresponding
#' to competition), else only one triangle will be set to a negative value (corresponding
#' to parasitism, predation or amensalism). If symmetric is true, the number of interactions
#' is counted only in one triangle of the interaction matrix.
#'
#' @param A interaction matrix
#' @param perc percentage of interactions to be modified
#' @param negval the negative value to be set
#' @param only introduce symmetric interactions
#' @return the modified interaction matrix
#' @export

setNeg<-function(A,perc=50, negval=-1, symmetric=TRUE){
  # edge number: number of non-zero entries in the interaction matrix
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
    A[xpos,ypos]=negval
    if(symmetric == TRUE){
      A[ypos,xpos]=negval
    }
    #print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
    #print(paste("reverse value:",A[ypos,xpos],sep=" "))
  }
  return(A)
}
