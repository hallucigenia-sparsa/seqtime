#' @title Plot mean variance versus the number of time points.
#'
#' @description Plot the mean variance across all taxa for an increasing number of time points. If groups are given,
#' compute the mean of beta diversity (assessed with Bray Curtis) across samples for each time point.
#' Groups are interpreted as replicates, such that each group is assumed to cover the same time points in the same order.
#'
#' @param x a community matrix with rows representing taxa and columns time points
#' @param groups a vector with group assignments and as many entries as there are samples
#' @param plot do the plot
#' @param sd use the standard deviation instead of the variance
#' @return a list with variances (or standard deviations if sd is true) and the intersection, slope, pval and adjR2 of the linear regression of variances against time
#' @examples
#' N=20
#' ts=glv(N=N,generateA(N))
#' out=varEvol(ts)
#' @export

varEvol<-function(x, groups=c(), plot=TRUE, sd=FALSE){
  meancumvar=c()
  var.measure="variance"
  if(sd==TRUE){
      var.measure="standard deviation"
  }
  # replicates provided
  if(length(groups)>0){
    var.measure="mean"
    group.ids=unique(groups)
    number.timepoints=length(which(groups==group.ids[1]))
    print(paste("Number of replicates:",length(group.ids)))
    if(length(group.ids)<2){
      stop("At least two groups required.")
    }
    #print(paste("Replicates are supposed to have",number.timepoints,"time points each."))
    values.per.group=list()
    for(group.id in group.ids){
      values.per.group[[group.id]]=x[,which(groups==group.id)]
      if(ncol(values.per.group[[group.id]]) != number.timepoints){
        stop(paste("Group",group.id,"does not have as many time points as the first group!"))
      }
    }
    # loop time points
    for(tp in 1:number.timepoints){
      replicate.mat=matrix(NA,nrow=nrow(x),ncol=length(group.ids)) # as many replicates as groups
      # collect samples at the same time point across replicates
      for(group.index in 1:length(group.ids)){
        replicate.mat[,group.index]=values.per.group[[group.ids[group.index]]][,tp]
      }
      # compute all pair-wise beta dissimilarities
      # vegdist computes dissimilarities row-wise, here should be computed column-wise
      bc.dist.mat=as.matrix(vegdist(t(replicate.mat)))
      #print(dim(bc.dist.mat))
      # get lower triangle elements from symmetric matrix (so upper triangle works equally well)
      bc.dist=bc.dist.mat[lower.tri(bc.dist.mat)]
      meancumvar=c(meancumvar,mean(bc.dist))
    } # end loop time points
  }else{
    cumvar=matrix(nrow=nrow(x),ncol=(ncol(x)-1))
    for(tp in 2:ncol(x)){
      for(i in 1:nrow(x)){
        if(sd==TRUE){
          cumvar[i,tp-1]=sd(x[i,1:tp])
        }else{
          cumvar[i,tp-1]=var(x[i,1:tp])
        }
      }
    }
    meancumvar=apply(cumvar,2,mean)
  }

  time=c(1:length(meancumvar))

  # do plot
  if(plot==TRUE){
    if(length(groups)==0){
      plot(meancumvar,xlab="Time",ylab="",main=paste("Change of ",var.measure," over time", sep=""))
      title(ylab=paste("Mean ",var.measure," up to time point", sep=""), line=4.5)
    }else{
      plot(meancumvar,xlab="Time",ylab=paste(firstup(var.measure)," of Bray Curtis dissimilarity",sep=""),main=paste("Change of ",var.measure," of dissimilarity across replicates", sep=""))
    }
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

# taken from: http://stackoverflow.com/questions/18509527/first-letter-to-upper-case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
