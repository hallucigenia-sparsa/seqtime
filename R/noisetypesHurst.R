#' Plot noise types versus their range of Hurst exponents
#'
#' Given a matrix and the output of noise type classification, compute the
#' Hurst exponent of each matrix row and draw a box plot that plots for each
#' noise type the range of Hurst exponents. Hurst exponents are computed with
#' function HurstK in package FGN.
#' @param x a matrix
#' @param noisetypes the noisetypes of the matrix
#' @param header header string
#' N=50
#' M=500
#' metapop=generateAbundances(N=M, mode=5, probabs=TRUE)
#' ts=simHubbell(N=N, M=M,I=1500,d=N, m.vector=metapop, tskip=50, tend=500)
#' noisetypes=identifyNoisetypes(ts)
#' noisetypesHurst(ts,noisetypes,header="Hurst exponent stratified by noise type")
#' @export

noisetypesHurst<-function(x, noisetypes, header=""){
  hurst.pink=c()
  if(length(noisetypes$pink) > 0){
    for(i in 1:length(noisetypes$pink)){
      h=FGN::HurstK(x[noisetypes$pink[i],])
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
      h=FGN::HurstK(x[noisetypes$brown[i],])
      hurst.brown=c(hurst.brown,h)
    }
  }
  hurst.black=c()
  if(length(noisetypes$black) > 0){
    for(i in 1:length(noisetypes$black)){
      h=FGN::HurstK(x[noisetypes$black[i],])
      hurst.black=c(hurst.black,h)
    }
  }
  print(paste("white mean Hurst:",mean(hurst.white)))
  print(paste("pink mean Hurst:",mean(hurst.pink)))
  print(paste("brown mean Hurst:",mean(hurst.brown)))
  print(paste("black mean Hurst:",mean(hurst.black)))

  hursts=list(hurst.white,hurst.pink,hurst.brown,hurst.black)
  names(hursts)=c("white","pink","brown","black")
  boxplot(hursts,col=c("white","pink","brown","black"),ylab="Hurst exponent", main=header)
}
