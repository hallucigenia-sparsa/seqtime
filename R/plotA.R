#' @title Plot an interaction matrix.
#'
#' @description By default, the interaction strengths are set to -1 or 1. Negative values
#' are plotted in red, positive in green. If scale.weight is true, the interaction
#' strengths are scaled to lie in the range of [-1,1]. The option original suppresses
#' any modification of the interaction strengths. If interaction strengths are scaled or the original
#' ones are used, the method ggplot is recommended, since it adds a color legend.
#' Alternatively, the interaction matrix can also be visualized as a directed network with method network.
#' In this case, the igraph network object is returned, to allow for manual adjustment of the layout with igraph's tkplot function.
#' Method ggplot requires ggplot2 and reshape2.
#'
#' @param A interaction matrix
#' @param method image, ggplot or network (ggplot requires ggplot2 and reshape2, image is therefore default), image and ggplot will plot the matrix, network will plot the network with igraph
#' @param header the title of the plot (does not apply to method network)
#' @param scale.weight scale interaction strengths between -1 and 1 (does not apply to method network)
#' @param original plot original values (does not apply to method network)
#' @param setNAToZero set missing values to zeros
#' @param removeOrphans remove orphan nodes (method network)
#' @param removeLoops remove loops (method network)
#' @param returnNetwork return the network for manual adjustment with tkplot (method network)
#' @examples
#' plotA(generateA(20,c=0.1))
#' @export

plotA<-function(A, method="image", header="", scale.weight=FALSE, original=FALSE, setNAToZero=FALSE, removeOrphans=TRUE, removeLoops=FALSE, returnNetwork=FALSE){
  A=as.matrix(A)
  if(setNAToZero==TRUE){
    A[is.na(A)]=0
  }
  if(method=="network"){
    scale.weight=FALSE
    original=TRUE
  }
  old.par=par()
  # scale values from -1 to 1
  max=max(A,na.rm=TRUE)
  min=min(A,na.rm=TRUE)
  print(paste("Largest value:",max))
  print(paste("Smallest value:",min))
  min=abs(min)
  if(original == FALSE){
    for(i in 1:nrow(A)){
      for(j in 1:ncol(A)){
        val=A[i,j]
        if(!is.na(val) && val > 0){
          if(scale.weight == TRUE){
            val=val/max
          }else{
            val=1
          }
        }else if(!is.na(val) && val < 0){
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
    # adjust palette, else in the absence of positive edges, interactions will appear red on a green bg
    if(length(A[A>0])==0){
      palette <- colorRampPalette(c("red","white"))
      colorNumber=2
    }
    if(length(A[A<0])==0){
      palette <- colorRampPalette(c("white","green"))
      colorNumber=2
    }
    if(scale.weight == TRUE){
      colorNumber=40
    }
    image(A,col=palette(colorNumber),main=header,axes=TRUE,xaxs="r",yaxs="r")
  }else if(method == "ggplot"){
    # https://stat.ethz.ch/pipermail/r-help/2005-September/078958.html
    searchGgplot2=length(grep(paste("^package:","ggplot2", "$", sep=""), search()))
    searchReshape2=length(grep(paste("^package:","reshape2", "$", sep=""), search()))
    # check whether ggplot2 is there
    if (searchGgplot2==0) {
      stop("ggplot2 is not installed/loaded. It is required for this option.")
    }
    # check whether reshape2 is there
    if (searchReshape2==0) {
      stop("reshape2 is not installed/loaded. It is required for this option.")
    }
    scale.plot<-max(c(max(A),-min(A)))
    Var1=NA
    Var2=NA
    value=NA
    # theme(axis.text.x = element_text(angle = 90, hjust = 1))+ coord_fixed() + coord_fixed()
    p1<-ggplot2::ggplot(reshape2::melt(A), ggplot2::aes(Var1,Var2, fill=value)) + ggplot2::geom_raster()+ ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "green",limits=c(-scale.plot, scale.plot)) + ggplot2::ggtitle(header) + ggplot2::labs(x = "",y="")
    plot(p1)
  }else if(method=="network"){
    if(removeOrphans==TRUE){
      A=modifyA(A,mode="removeorphans")
    }
    g=graph.adjacency(A, mode="directed",weighted=TRUE)
    if(removeLoops==TRUE){
      g=simplify(g, remove.multiple = F, remove.loops = T)
    }
    colors=c()
    for(weight in E(g)$weight){
      if(is.na(weight)){
        colors=append(colors,"gray")
      }else if(weight<0){
        colors=append(colors,"red")
      }else if(weight>0){
        colors=append(colors,"green")
      }
    }
    # attributes can be placed in plot also, but this has the advantage of easier transfer to tkplot, see http://kateto.net/network-visualization
    E(g)$arrow.size=0.3
    E(g)$color=colors
    V(g)$color="white"
    V(g)$frame.color="black"
    V(g)$label.color="black"
    # alternative: https://www.ggplot2-exts.org/ggraph.html
    # http://kateto.net/network-visualization
    plot(g,layout=layout.fruchterman.reingold(g))
    if(returnNetwork==TRUE){
      return(g)
    }
  }
  par=old.par
}
