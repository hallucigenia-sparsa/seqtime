#' Plot an interaction matrix.
#'
#' For plotting, the interaction strengths are scaled to lie between -1 and 1. Negative values
#' are plotted in red, positive in green.
#' If corrplot is available, it can be used for plotting.
#'
#' @param A interaction matrix
#' @param method image or corrplot (corrplot requires the corrplot package being loaded first, image is therefore default)
#' @param header the title of the plot
#' @param display.weight color-code the interaction strength by varying the red/green hues (only for corr.plot)
#' @param show.grid display the grid (only for corrplot)
#' @examples
#' plotA(generateA(20,c=0.1))
#' @export

plotA<-function(A, method="image",header="",display.weight=FALSE, show.grid=FALSE){
  A=as.matrix(A)

  # scale values from -1 to 1 and use corrplot
  # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  max=max(A)
  min=min(A)
  print(paste("Largest value:",max))
  print(paste("Smallest value:",min))
  min=abs(min)
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      val=A[i,j]
      if(val > 0){
        if(display.weight == TRUE){
          val=val/max
        }else{
          val=1
        }
      }else if(val < 0){
        if(display.weight == TRUE){
          val = val/min
        }else{
          val=-1
        }
      }
      A[i,j]=val
    } # column loop
  } # row loop

  if(method == "image"){
    palette <- colorRampPalette(c("red","white","green"))
    heatmap(A,Colv=NA,Rowv=NA,keep.dendro=FALSE,reorderfun=NULL,col=palette(3),main=header)
  }else if(method == "corrplot"){
    # check whether corrplot is there and load if it is installed
    # thanks to: http://r-pkgs.had.co.nz/description.html
    if (!requireNamespace("corrplot", quietly = TRUE)) {
      stop("corrplot is not installed. Please install it.", call. = FALSE)
    }
    grid.color="white"
    if(show.grid == TRUE){
      grid.color="gray"
    }
    if(display.weight){
      corrplot::corrplot(A,title=header,method="square",is.corr=TRUE,tl.pos = "n",bg="white",addgrid.col = grid.color)
    }else{
      # tl.pos=n: suppress text label
      corrplot::corrplot(A,title=header,method="square",col=c("red","white","green"),addgrid.col = grid.color,bg="white",is.corr=TRUE, tl.pos = "n")
    }
  }
}


