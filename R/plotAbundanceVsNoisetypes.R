#' @title Plot taxon abundances classified by noise type
#'
#' @description A box plot with taxon abundances is drawn for each noise type.
#'
#' @param x a matrix with taxa as rows and samples as columns
#' @param noisetypes the noisetypes object generated for the matrix
#' @export

plotAbundanceVsNoisetypes<-function(x, noisetypes){
  white.noise=c()
  pink.noise=c()
  brown.noise=c()
  black.noise=c()
  # pool counts of taxa with the same noise type
  if(length(noisetypes$white)>0){
    white.noise=as.vector(x[noisetypes$white,])
  }else{
    white.noise=c(NA)
  }
  if(length(noisetypes$pink)>0){
    pink.noise=as.vector(x[noisetypes$pink,])
  }else{
    pink.noise=c(NA)
  }
  if(length(noisetypes$brown)>0){
    brown.noise=as.vector(x[noisetypes$brown,])
  }else{
    brown.noise=c(NA)
  }
  if(length(noisetypes$black)>0){
    black.noise=as.vector(x[noisetypes$black,])
  }else{
    black.noise=c(NA)
  }
  noise.types.bars=list(white.noise,pink.noise,brown.noise,black.noise)
  names(noise.types.bars)=c("white noise","pink noise","brown noise","black noise")
  boxplot(noise.types.bars, ylab="Abundance",main="Abundance versus noise type", col=c("white","pink","brown","black"))
}
