#' Plot an interaction matrix.
#'
#' By default, the interaction strengths are set to -1 or 1. Negative values
#' are plotted in red, positive in green. If scale.weight is true, the interaction
#' strengths are scaled to lie in the range of [-1,1]. The option original suppresses
#' any modification of the interaction strengths. If interaction strengths are scaled or the original
#' ones are used, the method ggplot is recommended, since it adds a color legend.
#' Method ggplot requires ggplot2 and reshape2.
#'
#' @param A interaction matrix
#' @param method image or ggplot (ggplot requires ggplot2 and reshape2, image is therefore default)
#' @param header the title of the plot
#' @param scale.weight scale interaction strengths between -1 and 1
#' @param original plot original values
#' @examples
#' plotA(generateA(20,c=0.1))
#' @export

plotA<-function(A, method="image", header="", scale.weight=FALSE, original=FALSE){
  A=as.matrix(A)

  old.par=par()
  # scale values from -1 to 1
  max=max(A)
  min=min(A)
  print(paste("Largest value:",max))
  print(paste("Smallest value:",min))
  min=abs(min)
  if(original == FALSE){
    for(i in 1:nrow(A)){
      for(j in 1:ncol(A)){
        val=A[i,j]
        if(val > 0){
          if(scale.weight == TRUE){
            val=val/max
          }else{
            val=1
          }
        }else if(val < 0){
          if(scale.weight == TRUE){
            val = val/min
          }else{
            val=-1
          }
        }
        A[i,j]=val
      } # column loop
    } # row loop
  }
  if(method == "image"){
    palette <- colorRampPalette(c("red","white","green"))
    colorNumber=3
    if(scale.weight == TRUE){
      colorNumber=40
    }
    image(A,col=palette(colorNumber),main=header,axes=TRUE,xaxs="r",yaxs="r")
  }else if(method == "ggplot"){
    # check whether ggplot2 is there
    if (!require("ggplot2")) {
      stop("ggplot2 is not installed. Please install it.")
    }
    # check whether reshape2 is there
    if (!require("reshape2")) {
      stop("reshape2 is not installed. Please install it.")
    }
    scale.plot<-max(c(max(A),-min(A)))
    # theme(axis.text.x = element_text(angle = 90, hjust = 1))+ coord_fixed() + coord_fixed()
    p1<-ggplot2::ggplot(reshape2::melt(A), aes(Var1,Var2, fill=value)) + geom_raster()+ scale_fill_gradient2(low = "red", mid = "white", high = "green",limits=c(-scale.plot, scale.plot)) + ggtitle(header) + labs(x = "",y="")
    plot(p1)
  }
  par=old.par
}


