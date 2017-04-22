#' @title Do a barplot of the noise type counts
#'
#' @description Plot the counts for each noise type in a barplot.
#'
#' @param noise the noisetypes object
#' @rdname plot
#' @method plot noisetypes
#' @export
plot.noisetypes<-function(noise){
  noise.types.bars=c(length(noise$white),length(noise$pink),length(noise$brown), length(noise$black))
  names(noise.types.bars)=c("white noise","pink noise","brown noise", "black noise")
  barplot(noise.types.bars, col=c("white","pink","brown","black"),main="Noise types", ylab="Number of rows")
}
