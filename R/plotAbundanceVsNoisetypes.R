#' @title Plot taxon abundances classified by noise type
#'
#' @description A box plot with taxon abundances is drawn for each noise type.
#'
#' @param x a matrix with taxa as rows and samples as columns
#' @param nt the noisetypes object generated for the matrix
#' @param sig assess significance of difference between mean abundances with ANOVA and Tukey's post-hoc test
#' @export

plotAbundanceVsNoisetypes<-function(x, nt, sig=FALSE){
  white.noise=c()
  mean.abundance.white=c()
  pink.noise=c()
  mean.abundance.pink=c()
  brown.noise=c()
  mean.abundance.brown=c()
  black.noise=c()
  mean.abundance.black=c()
  # pool counts of taxa with the same noise type
  if(length(nt$white)>0){
    white.noise=as.vector(x[nt$white,])
    mean.abundance.white=rowMeans(x[nt$white,])
  }else{
    white.noise=c(NA)
  }
  if(length(nt$pink)>0){
    pink.noise=as.vector(x[nt$pink,])
    mean.abundance.pink=rowMeans(x[nt$pink,])
  }else{
    pink.noise=c(NA)
  }
  if(length(nt$brown)>0){
    brown.noise=as.vector(x[nt$brown,])
    mean.abundance.brown=rowMeans(x[nt$brown,])
  }else{
    brown.noise=c(NA)
  }
  if(length(nt$black)>0){
    black.noise=as.vector(x[nt$black,])
    mean.abundance.black=rowMeans(x[nt$black,])
  }else{
    black.noise=c(NA)
  }
  noise.types.bars=list(white.noise,pink.noise,brown.noise,black.noise)
  if(sig==TRUE){
    mean.abundances=c(mean.abundance.white,mean.abundance.pink,mean.abundance.brown,mean.abundance.black)
    noisetype.class=c(rep("white",length(nt$white)),rep("pink",length(nt$pink)),rep("brown",length(nt$brown)),rep("black",length(nt$black)))
    nt.data=data.frame(mean.abundances,noisetype.class)
    nt.data$noisetype.class=as.factor(nt.data$noisetype.class)
    nt.lm=lm(mean.abundances ~ noisetype.class, data=nt.data)
    nt.av=aov(nt.lm)
    tukey.test=TukeyHSD(nt.av)
    print(tukey.test)
    #print(wilcox.test(mean.abundance.brown,mean.abundance.pink))
  }
  names(noise.types.bars)=c("white noise","pink noise","brown noise","black noise")
  boxplot(noise.types.bars, ylab="Abundance",main="Abundance versus noise type", col=c("white","pink","brown","black"))
}
