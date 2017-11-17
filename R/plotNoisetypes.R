#' @title Do a barplot of the noise types
#'
#' @description Plot the count for each noise type in a barplot.
#'
#' @param noisetypes the noisetypes object
#' @param percentage plot the percentage instead of the number
#' @param total total number of items (in case not all items have a noise type assigned, needed to compute percentages with reference to all items)
#' @param stacked plot a stacked barplot
#' @param sample.id the sample identifier (optional)
#' @export
plotNoisetypes<-function(noisetypes, percentage=TRUE, total=NA, stacked=FALSE, sample.id=""){
  noise.types.bars=c(length(noisetypes$white),length(noisetypes$pink),length(noisetypes$brown), length(noisetypes$black))
  if(percentage){
    if(is.na(total)){
      total=sum(noise.types.bars)
    }
    onePerc=total/100
    noise.types.bars=noise.types.bars/onePerc
  }
  if(stacked){
    noise.types.bars=cbind(noise.types.bars,rep(NA,length(noise.types.bars)))
    colnames(noise.types.bars)=c(sample.id,"")
  }else{
    names(noise.types.bars)=c("white noise","pink noise","brown noise", "black noise")
  }
  ylab="Number of taxa"
  if(percentage){
    ylab="Percentage of taxa"
  }
  main="Noise types"
  if(sample.id!=""){
    main=paste(main,sample.id)
  }
  barplot(noise.types.bars, col=c("white","pink","brown","black"),main=main, ylab=ylab)
}
