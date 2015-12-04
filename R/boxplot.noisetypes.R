#' Do a boxplot of matrix rows classified by noise type
#'
#' @param x a matrix
#' @param noisetypes the noisetypes object generated for the matrix
#' @export

boxplot.noisetypes<-function(x, noisetypes){
  # pool counts of taxa with the same noise type
  white.noise=as.vector(x[noisetypes$white,])
  pink.noise=as.vector(x[noisetypes$pink,])
  brown.noise=as.vector(x[noisetypes$brown,])
  noise.types.bars=list(white.noise,pink.noise,brown.noise)
  names(noise.types.bars)=c("white noise","pink noise","brown noise")
  boxplot(noise.types.bars)
}
