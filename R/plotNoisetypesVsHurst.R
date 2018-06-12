#' @title Plot noise types versus their range of Hurst exponents
#'
#' @description Given a matrix and the output of noise type classification, compute the
#' Hurst exponent of each matrix row and draw a box plot that plots for each
#' noise type the range of Hurst exponents. Hurst exponents are computed with
#' function hurstexp in the pracma package (simple R/S Hurst is used).
#' @param x a matrix
#' @param noisetypes the noisetypes of the matrix
#' @param d window size for Hurst exponent computation with function hurstexp
#' @param header header string
#' @examples
#' \dontrun{
#' N=50
#' M=500
#' metapop=generateAbundances(N=M, mode=5, probabs=TRUE)
#' ts=simHubbell(N=N, M=M,I=1500,d=N, m.vector=metapop, tskip=50, tend=500)
#' noisetypes=identifyNoisetypes(ts)
#' noisetypesHurst(ts,noisetypes,header="Hurst exponent stratified by noise type")
#' }
#' @export

plotNoisetypesVsHurst<-function(x, noisetypes,d=round(ncol(x)/3), header=""){
  hurst.pink=c()
  if(length(noisetypes$pink) > 0){
    for(i in 1:length(noisetypes$pink)){
      h=pracma::hurstexp(x[noisetypes$pink[i],],d=d)$Hs
      hurst.pink=c(hurst.pink,h)
    }
  }else{
    hurst.pink=c(NA)
  }
  hurst.white=c()
  if(length(noisetypes$white) > 0){
    for(i in 1:length(noisetypes$white)){
      h=pracma::hurstexp(x[noisetypes$white[i],],d=d)$Hs
      hurst.white=c(hurst.white,h)
    }
  }else{
    hurst.white=c(NA)
  }
  hurst.brown=c()
  if(length(noisetypes$brown) > 0){
    for(i in 1:length(noisetypes$brown)){
      h=pracma::hurstexp(x[noisetypes$brown[i],],d=d)$Hs
      hurst.brown=c(hurst.brown,h)
    }
  }else{
    hurst.brown=c(NA)
  }
  hurst.black=c()
  if(length(noisetypes$black) > 0){
    for(i in 1:length(noisetypes$black)){
      h=pracma::hurstexp(x[noisetypes$black[i],],d=d)$Hs
      hurst.black=c(hurst.black,h)
    }
  }else{
    hurst.black=c(NA)
  }
  print(paste("white mean Hurst:",mean(hurst.white)))
  print(paste("pink mean Hurst:",mean(hurst.pink)))
  print(paste("brown mean Hurst:",mean(hurst.brown)))
  print(paste("black mean Hurst:",mean(hurst.black)))

  hursts=list(hurst.white,hurst.pink,hurst.brown,hurst.black)
  names(hursts)=c("white","pink","brown","black")
  boxplot(hursts,col=c("white","pink","brown","black"),ylab="Hurst exponent", main=header)
}
