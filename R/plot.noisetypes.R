#' @title Do a barplot of the noise type counts
#'
#' @param noisetypes the noisetypes object
#' @export

plot.noisetypes<-function(noisetypes){
  noise.types.bars=c(length(noisetypes$white),length(noisetypes$pink),length(noisetypes$brown), length(noisetypes$black))
  names(noise.types.bars)=c("white noise","pink noise","brown noise", "black noise")
  barplot(noise.types.bars, col=c("white","pink","brown","black"),main="Noise types", ylab="Number of rows")
}
