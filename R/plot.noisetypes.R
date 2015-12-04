#' Do a barplot of the noise type counts
#'
#' @param noisetypes the noisetypes object
#' @export

plot.noisetypes<-function(noisetypes){
	noise.types.bars=c(length(noisetypes$white),length(noisetypes$pink),length(noisetypes$brown))
	names(noise.types.bars)=c("white noise","pink noise","brown noise")
	barplot(noise.types.bars, col=c("white","pink","brown"),main="Noise types", ylab="Number of rows")
}