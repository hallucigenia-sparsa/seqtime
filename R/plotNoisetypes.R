#' @title Do a barplot of the noise types
#'
#' @description Plot the count for each noise type in a barplot.
#'
#' @param noisetypes the noisetypes object
#' @export
plotNoisetypes<-function(noisetypes){
  noise.types.bars=c(length(noisetypes$white),length(noisetypes$pink),length(noisetypes$brown), length(noisetypes$black))
  names(noise.types.bars)=c("white noise","pink noise","brown noise", "black noise")
  barplot(noise.types.bars, col=c("white","pink","brown","black"),main="Noise types", ylab="Number of rows")
}
