#' Plot noise types versus their range of Hurst exponents
#'
#' Given a matrix and the output of noise type classification, compute the
#' Hurst exponent of each matrix row and draw a box plot that plots for each
#' noise type the range of Hurst exponents. Hurst exponents are computed with
#' function HurstK in package FGN.
#' @param x a matrix
#' @param noisetypes the noisetypes of the matrix
#' @param header header string
#' @export

noisetypesHurst<-function(x, noisetypes, header=""){
  hurst.pink=c()
  if(length(noisetypes$pink) > 0){
    for(i in 1:length(noisetypes$pink)){
      h=HurstK(x[noisetypes$pink[i],])
      hurst.pink=c(hurst.pink,h)
    }
  }
  hurst.white=c()
  if(length(noisetypes$white) > 0){
    for(i in 1:length(noisetypes$white)){
      h=FGN::HurstK(x[noisetypes$white[i],])
      hurst.white=c(hurst.white,h)
    }
  }
  hurst.brown=c()
  if(length(noisetypes$brown) > 0){
    for(i in 1:length(noisetypes$brown)){
      h=HurstK(x[noisetypes$brown[i],])
      hurst.brown=c(hurst.brown,h)
    }
  }
  print(paste("pink mean Hurst:",mean(hurst.pink)))
  print(paste("brown mean Hurst:",mean(hurst.brown)))
  print(paste("white mean Hurst:",mean(hurst.white)))

  hursts=list(hurst.white,hurst.brown,hurst.pink)
  names(hursts)=c("white", "brown", "pink")
  boxplot(hursts,col=c("white","brown","pink"),ylab="Hurst exponent", main=header)
}
