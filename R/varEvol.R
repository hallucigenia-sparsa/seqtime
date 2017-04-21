#' @title Plot change of mean variance versus the number of time points.
#'
#' @param x a community matrix with rows representing taxa and columns time points
#' @param plot do the plot
#' @return a list with variances and the intersection, slope, pval and adjR2 of the linear regression of variances against time
#' @examples
#' N=20
#' ts=glv(N=N,generateA(N))
#' out=varEvol(ts)
#' @export

varEvol<-function(x, plot=TRUE){
  cumvar=matrix(nrow=nrow(x),ncol=(ncol(x)-1))
  for(tp in 2:ncol(x)){
    for(i in 1:nrow(x)){
      cumvar[i,tp-1]=var(x[i,1:tp])
    }
  }
  meancumvar=apply(cumvar,2,mean)
  time=c(1:length(meancumvar))

  # do plot
  if(plot==TRUE){
    plot(meancumvar,xlab="Time",ylab="",main="Change of variance over time")
    title(ylab="Mean variance up to time point", line=4.5)
  }

  # do linear regression
  linreg = lm(formula = meancumvar~time)
  intersection = linreg$coefficients[1]
  slope=linreg$coefficients[2]
  sum=summary(linreg)
  adjR2=sum$adj.r.squared
  pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
  res=list(meancumvar, intersection, slope, pval, adjR2)
  names(res)=c("variances", "intersection","slope","pval","adjR2")
  return(res)
}
