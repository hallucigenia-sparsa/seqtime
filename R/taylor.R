#' @title Plot relationship between row mean and row variance
#'
#' @description The power law relationship between variance and mean is known as Taylor's law, which is
#' defined as: \eqn{var(Y) = a*mean(Y)^b}. We can obtain the power as the slope in log-scale: \eqn{log(var(Y)) = log(a)+b*log(mean(Y))}
#' @details Confidence intervals are estimated by first resampling columns with replacement from the input matrix and then recomputing log means
#' and log variances on the bootstrapped data.
#'
#' @param x a matrix, mean and variance are computed row-wise
#' @param type the type of plot to do: mean.var (mean vs variance), boxplot (row-wise), taylor (powerlaw fitted to mean vs variance)
#' @param boot compute confidence interval through bootstrapping with given number of iterations (only for type=taylor)
#' @param lower.conf lower confidence interval
#' @param upper.conf upper confidence interval
#' @param pseudo add a pseudo count to deal with zeros in log-log plot (for type=taylor)
#' @param col the color of the dots
#' @param header header string
#' @param label label points interactively (only for type=taylor, escape to stop labeling)
#' @param plot whether to do the plot
#' @return for type taylor, the slope, p-value, adjusted R2 of the Taylor law as well as the log means and variances are returned (slope, pval, adjR2, logmeans, logvars)
#' @references L.R. Taylor (1961). Aggregation, variance and the mean. Nature 189, 732-735.
#' @export

taylor<-function(x, type="taylor", boot=1000, lower.conf=0.05, upper.conf=0.95, pseudo=0, col="black", header="", label=FALSE, plot=TRUE){
  pval=NA
  slope=NA
  adjR2=NA

  # compute means and variances
  if(type=="taylor" || type=="mean.var"){
      means=apply(x,1, mean, na.rm=TRUE)
      vars=apply(x,1, var, na.rm=TRUE)
      logvars=log(vars+pseudo)
      logmeans=log(means+pseudo)
      #print(paste("ori logvars",length(logvars)))
      lowerConfV=c()
      upperConfV=c()
      lowerConfM=c()
      upperConfM=c()
      # bootstraps are only carried out for Taylor law
      if(!is.na(boot) && boot>0 && type=="taylor"){
        #bootstrapped.line=matrix(NA, nrow=boot, ncol=nrow(x)) # iterations x rows
        bootstrapped.means=matrix(NA, nrow=boot, ncol=nrow(x)) # iterations x rows
        bootstrapped.vars=matrix(NA, nrow=boot, ncol=nrow(x)) # iterations x rows
        for(iteration in 1:boot){
          # generate bootstrapped data
          x.b.indices=sample(c(1:ncol(x)),size=ncol(x),replace = TRUE)
          x.b=x[,x.b.indices]
          # compute Taylor law for bootstrapped data
          b.iter.means=apply(x.b,1, mean, na.rm=TRUE)
          b.iter.vars=apply(x.b,1, var, na.rm=TRUE)
          b.iter.logvars=log(b.iter.vars+pseudo)
          b.iter.logmeans=log(b.iter.means+pseudo)
          bootstrapped.means[iteration,]=b.iter.logmeans
          bootstrapped.vars[iteration,]=b.iter.logvars
          #b.linreg = lm(formula = b.iter.logvars~b.iter.logmeans)
          #outL=lowess(b.iter.logmeans,b.iter.logvars)
          #print(paste("fitted logvars",length(b.linreg$fitted.values)))
          #bootstrapped.line[iteration,]=b.linreg$fitted.values # keep fitted taylor law
        } # end bootstrap iterations
        # get confidence intervals for each row
        for(row.index in 1:nrow(x)){
          lowerConfV=c(lowerConfV, quantile(bootstrapped.vars[,row.index],lower.conf))
          upperConfV=c(upperConfV, quantile(bootstrapped.vars[,row.index],upper.conf))
          lowerConfM=c(lowerConfM, quantile(bootstrapped.means[,row.index],lower.conf))
          upperConfM=c(upperConfM, quantile(bootstrapped.means[,row.index],upper.conf))
        }
      } # end do bootstrap
  }

  # plot
  if(type == "boxplot"){
    #par(las=3, cex=0.5)
    if(plot == TRUE){
      boxplot(t(x), ylab="Abundances", xaxt='n', ann=FALSE, border=col)
    }
  }else{
    if(type == "taylor"){
      reg.data=data.frame(logvars,logmeans)
      linreg = lm(formula = logvars~logmeans)
      # print(paste("Intercept:",linreg$coefficients[1]))
      # print(paste("Slope:",linreg$coefficients[2]))
      slope=linreg$coefficients[2]
      sum=summary(linreg)
      adjR2=sum$adj.r.squared
      #print(paste("Adjusted R2:",sum$adj.r.squared))
      pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
      # print(paste("P-value:",pval))
      if(plot == TRUE){
        xlim=c(min(logmeans), max(logmeans))  # -1, +2
        ylim=c(min(logvars), max(logvars))
        plot(logmeans,logvars, xlim=xlim, ylim=ylim,xlab="Log(mean)", ylab="Log(variance)", main=paste("Taylor's law",header,"\np-value:",round(pval,3),", adjusted R2:",round(sum$adj.r.squared,2),", slope:",round(linreg$coefficients[2],2)),type="p",pch="+",col=col)
        abline(linreg,bty="n",col="red")
        if(!is.na(boot) && boot>0){
          fillcol=rgb(1,0,0,alpha=0.5)
          # error bars in both directions:
          for(val.index in 1:length(logvars)){
            arrows(logmeans[val.index],upperConfV[val.index], logmeans[val.index],lowerConfV[val.index], col=fillcol, length=0.05, angle=90)
            arrows(logmeans[val.index],lowerConfV[val.index], logmeans[val.index],upperConfV[val.index], col=fillcol, length=0.05, angle=90)
            arrows(upperConfM[val.index],logvars[val.index], lowerConfM[val.index],logvars[val.index], col=fillcol, length=0.05, angle=90)
            arrows(lowerConfM[val.index],logvars[val.index], upperConfM[val.index],logvars[val.index], col=fillcol, length=0.05, angle=90)
          }
        }
        if(label){
          #text(logmeans,logvars,labels=rownames(x), pos=3, cex=0.9)
          identify(logmeans,logvars,labels=rownames(x))
        }
      } # end plot true
      res=list(slope,pval,adjR2,logmeans,logvars)
      names(res)=c("slope","pval","adjR2","logmeans","logvars")
      return(res)
    }else if(type == "mean.var"){
      if(plot==TRUE){
        plot(means,vars,xlab="Mean", ylab="Variance", main="Mean versus variance",type="p",pch="+",col=col)
      }
      res=list(means,vars)
      names(res)=c("means","vars")
      return(res)
    } # mean variance end
  } # no box plot
}
