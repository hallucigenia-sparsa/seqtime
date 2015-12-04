#' Plot the time series row-wise
#'
#' @param x the matrix of time series
#' @param time.given if true, then the column names are supposed to hold the time units
#' @param num the number of rows to plot (starting from the first row)
#' @param legend add a legend
#' @param pch plot parameter
#' @param type pot parameter
#' @param header header string
#' @export

tsplot<-function(x, time.given=FALSE, num=nrow(x), legend=FALSE, pch="+", type="l", header=""){
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
  plot(time,as.numeric(x[1,]),ylim = range(x, na.rm = T),xlab = xlab, ylab = ylab, main = main, type = type, col = my.colors[1], pch=pch)
  # loop over rows in data
  for(i in 2:num){
    lines(time,as.numeric(x[i,]), col = my.colors[i], type = type, pch=pch)
  }
  if(legend == TRUE){
    legend("right",legend=rownames(x), lty = rep(1,nrow(x)), col = my.colors, merge = TRUE, bg = "white", text.col="black")
  }
}
