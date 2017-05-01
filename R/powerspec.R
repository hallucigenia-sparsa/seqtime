#' @title Report slope of periodogram in log scale
#'
#' @description The periodogram plots frequencies (f) versus their power (spectrum).
#' In case their relationship is well described by a line in log scale,
#' its slope can be used to determine the noise type of a time series.
#' If the slope is around -1, the time series displays 1/f (pink) noise.
#' If it is around -2, the time series displays 1/f^2 (brown) noise. If the slope
#' is even steeper, the time series displays black noise.
#'
#' @details The function uses stats::spectrum to compute the periodogram. It also reports the significance
#' and goodness of fit of the power law.
#'
#' @param v time series vector
#' @param plot plot the periodogram with the power law in log-scale
#' @param header header string
#' @param col color used in periodogram if plot is true
#' @return return the slope, p-value, adjusted R2, log frequencies and log spectra
#' @examples
#' brownNoise=cumsum(rnorm(500,mean=10))
#' out.spec=powerspec(brownNoise, header="brown noise", plot=TRUE)
#' @export

powerspec<-function(v, plot=FALSE, header="", col="blue"){
  # generates periodogram
  out=stats::spectrum(v, plot=FALSE)
  logfreq=log(out$freq)
  logspec=log(out$spec)
  reg.data=data.frame(logfreq,logspec)
  linreg = lm(formula = logspec~logfreq)
  intercept=linreg$coefficients[1]
  slope=linreg$coefficients[2]
  sum=summary(linreg)
  r2.adj=sum$adj.r.squared
  pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
  #print(paste("Slope:",slope))
  #print(paste("adjusted R2:",r2.adj))
  #print(paste("P-value:",pval))
  if(plot == TRUE){
    plot(logfreq,logspec,xlab="Log(Frequency)", ylab="Log(Spectrum)", main=paste("Periodogram",header,"\np-value:",round(pval,3),", adjusted R2:",round(sum$adj.r.squared,2),", slope:",round(linreg$coefficients[2],2)),type="p",pch="+",col=col)
    abline(linreg,bty="n",col="red")
  }
  res=list(slope,pval,r2.adj,logfreq,logspec)
  names(res)=c("slope","pval","adjR2","logfreq","logspec")
  return(res)
}
