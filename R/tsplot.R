#' @title Time Series Plot
#' @description Plot the time series row-wise.
#' @param x the matrix of time series
#' @param time.given if true, then the column names are supposed to hold the time units
#' @param num the number of rows to plot (starting from the first row)
#' @param sample.points sample.points
#' @param legend add a legend
#' @param header string added to the plot title
#' @param pch plot parameter (point character)
#' @param lty plot parameter (line type)
#' @param type plot parameter (plot type, e.g. dots, lines etc.)
#' @param \\dots Additional arguments passed to plot()
#' @export
tsplot <- function(x, time.given=FALSE, num=nrow(x), sample.points=c(), legend=FALSE, header="", pch="+", lty=1, type="l", ...){
  col.vec = seq(0,1,1/nrow(x))
  my.colors = hsv(col.vec)
  main=paste("Community time series",header)
  xlab="Time points"
  ylab="Abundance"
  if(time.given == TRUE){
    time=as.numeric(colnames(x))
  }else{
    time=c(1:ncol(x))
  }
  plot(time,as.numeric(x[1,]),ylim = range(x, na.rm = T),xlab = xlab, ylab = ylab, main = main, type = type, col = my.colors[1], pch=pch, lty=lty, ...)
  # loop over rows in data
  for(i in 2:num){
    lines(time,as.numeric(x[i,]), col = my.colors[i], type = type, pch=pch, lty=lty, ...)
  }
  # if non-empty, loop over sample points
  if(length(sample.points)>0){
    for(i in 1:length(sample.points)){
        abline(v=sample.points[i],col="gray", lty=lty, ...)
    }
  }
  if(legend == TRUE){
    legend("right",legend=rownames(x), lty = rep(1,nrow(x)), col = my.colors, merge = TRUE, bg = "white", text.col="black")
  }
}
