# Number of asymmetric negative interactions
#
# Get the number of asymmetric negative interactions
# in an interaction matrix.
# If amensalism is true, also those interactions
# count as asymmetric that are (-,0), else only
# (+,-) is counted.
# A interaction matrix
# amensalism include amensalistic interactions

getNumAsymNeg<-function(A, amensalism=TRUE){
  counter.asym.neg=0
  # first row loop
  for(i in 1:nrow(A)){
    # second row loop
    for(j in 2:i){
      if(amensalism == TRUE){
        if((A[i,j]>=0 && A[j,i]<0) || (A[j,i]>=0 && A[i,j]<0)){
          counter.asym.neg=counter.asym.neg + 1
        }
      }else{
        if((A[i,j]>0 && A[j,i]<0) || (A[j,i]>0 && A[i,j]<0)){
          counter.asym.neg=counter.asym.neg + 1
        }
      } # amensalism yes/no
    } # second row loop
  } # first row loop
  return(counter.asym.neg)
}
