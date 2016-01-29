#' Print interaction type statistics
#'
#' @param A the interaction matrix
#' @export

getAStats<-function(A){

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

  res=list(EX,P,Pm,Pc,Pe,Pp,Pa)
  names(res)=c("meanstrength","nbinteractions","nbmut","nbcomp","nbexp","nbcom","nbam")
  return(res)
}
