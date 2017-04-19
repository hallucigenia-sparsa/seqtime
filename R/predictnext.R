# Given current time step, predict taxon abundances at the next time step.
#
# If the community time series matrix x is given, predict the taxon abundances after lag time steps
# for all (t,t+lag) sample pairs. If memory m is larger than one, then m preceding samples are taken
# into account in a multiple regression: \eqn{x[,t+lag] = a + b1*x[,t] + b2*x[,t-1] + b3*[,t-2] + ... + bm*x[,t-m]}
# The output is a plot with the predictor time step on the x-axis and the adjusted R2 on the y-axis.
# Optionally, a time series vector y can be provided. In this case, y is predicted from the community
# time series as follows:
# \deqn{y(t+lag) = a + \sum\limits_{i=1}^n(bi*xi(t))}
# where i is the taxon index. Predictor taxa xi are selected with randomForest, thus the randomForest package is required.
# If y is provided, instead of the AIC or R2, the increase in prediction error is plotted for each taxon. The larger
# the increase (upon taxon permutation), the more important the taxon for correct prediction of y.
#
# param x a taxon matrix with rows representing taxa and columns samples
# param y optionally, a time series of values measured for the samples in x
# param time vector (needed for plotting, ignored when y is provided)
# param group.membership the group membership vector provides for each sample its group, samples in each group are supposed to form time series (only if y is provided)
# param lag how many time steps predictor abundances precede predicted abundances
# param method prediction method (linreg = linear regression with lm, ignored when y is provided)
# param fixstep keep the predictor at the first time step and only increase predicted step-wise, including lag (ignored when y is provided)
# param m take m previous samples into account (ignored when y is provided)
# param aic plot the AIC instead of the adjusted R2 (ignored when y is provided)
# param ntree the number of trees to use in random forest regression (only if y is provided)
# param plot switch off plotting
# param header a string that is attached to the plot title
# return If y is provided, the indices of all taxa with positive %incMSE (selected) as well as the residual sum of squares (rss) are returned.

predictnext<-function(x, y, time=c(1:ncol(x)), group.membership=c(), lag=1, method="linreg", fixstep=FALSE, m=1, aic=FALSE, ntree=5000, plot=TRUE, header=""){
  if(length(time) != ncol(x)){
    stop("The time vector needs as many entries as x has columns!")
  }
  adjr2s=c()
  aics=c()
  plottime=c()

  if(!is.null(y)){
    if(length(y) != ncol(x)){
      stop("Vector y does not have as many entries as x has columns!")
    }
    if (!require("randomForest")) {
      stop("randomForest is not installed. Please install it.")
    }
    indices.na=which(is.na(y))
    if(length(indices.na)>0){
      warning(paste("Vector y should not contain any missing values!",length(indices.na),"missing values are removed."))
      keep.indices=setdiff(c(1:length(y)),indices.na)
      y=y[keep.indices]
      x=x[,keep.indices]
      group.membership=group.membership[keep.indices]
    }
    # reformat the data such that group membership is taken into account
    if(length(groups)>0){
      groups=unique(group.membership)
      #print(paste("Groups:",groups))
      totalx=c()
      totaly=c()
      # loop over groups
      for(group in groups){
        member.indices=which(group.membership==group)
        if(length(member.indices)>lag){
          #print(paste("Group members",colnames(x)[member.indices]))
          # take only group members into account
          xgroup=x[,member.indices]
          totaly=c(totaly,y[member.indices[(1+lag):length(member.indices)]])
          # take the requested lag into account
          xgroup=x[,1:ncol(xgroup)-lag]
          #print(dim(xgroup))
          totalx=cbind(totalx,xgroup)
        }else{
          print(paste("Skipping group",group))
        }
      }
      # rows are taxa, columns samples, but transpose is done below
      print(paste("Number of taxa:",nrow(totalx)))
      print(paste("Number of samples combined:",ncol(totalx)))
      print(paste("Number of y values combined:",length(totaly)))
      rownames(totalx)=rownames(x)
      colnames(totalx)=as.character(c(1:ncol(totalx)))
      data=rbind(totaly,totalx)
    }else{
      if(lag == 0){
        data=rbind(y,x)
      }else{
        yl=y[(1+lag):length(y)]
        xl=x[,1:(ncol(x)-lag)]
        # avoid error because sample names do not match
        if(!is.vector(yl)){
          colnames(yl)=c(colnames(xl))
        }
        data=rbind(yl,xl)
      }
    }
    rownames(data)[1]="y"
    df=data.frame(t(data))
    y=df[["y"]]
    rf.out <- randomForest(df$y~.,data=df,importance=TRUE, proximity=TRUE, ntree=ntree)
    imp=importance(rf.out)
    # http://stats.stackexchange.com/questions/161709/interpreting-var-explained-in-random-forest-output
    # tested with ozone data set
    varExplained=rf.out$rsq[length(rf.out$rsq)]*100
    incMSE=imp[,1]
    # indices of taxa with positive incMSE
    selectedTaxa=which(incMSE>0)
    # prepare barplot of taxon importance for prediction of vector
    incMSE=na.omit(incMSE)
    selected=incMSE[incMSE>0]
    # residual sum of squares
    rss=sum((y-rf.out$predicted)^2)
    if(plot == TRUE){
      prevpar=par()
      par(las = 2, cex=0.8, mar = c(5, 20, 4, 2)) # rotate labels, decrease font size, enlarge margin
      barplot(selected,horiz=TRUE,xlab="%incMSE",main=paste("Predictors with positive %incMSE ",header,"\nRSS=",round(rss,2),", Variance explained=",varExplained,"%",sep=""))
      par=prevpar
    }
    # prepare result object
    result=list(selectedTaxa, rss)
    names(result)=c("selected","rss")
    return(result)
  }
  else{
    # loop over time steps
    for(step in 1:(ncol(x)-lag)){
      #print(paste("step",step))
      if(fixstep == TRUE){
        predictor = as.numeric(x[,1])
        predicted=as.numeric(x[,(step+lag)])
        linreg = lm(formula = predicted~predictor)
      }else{
        if(m > 1){
          if(step >= m){
            Y=as.numeric(x[,(step+lag)])
            predictor=matrix(nrow=nrow(x),ncol=m)
            colnames=c()
            formula="Y~"
            metadata.index=1
            # assemble predictor
            for(predict.index in (step-m+1):step){
              #print(paste("predict index",predict.index))
              #print(paste("access index",predict.index))
              predictor[,metadata.index]=as.numeric(x[,predict.index])
              colname=paste("X",metadata.index,sep="")
              if(predict.index == step){
                formula=paste(formula,colname,"",sep="")
              }else{
                formula=paste(formula,colname,"+",sep="")
              }
              colnames=c(colnames,colname)
              metadata.index=metadata.index+1
            }
            #print(paste("formula",formula))
            colnames(predictor)=colnames
            data=data.frame(Y,predictor)
            linreg = lm(formula, data)
          }else{
            # less than m samples, do nothing
          }
        }else{
          predictor=as.numeric(x[,step])
          predicted=as.numeric(x[,(step+lag)])
          linreg = lm(formula = predicted~predictor)
        }
      }
      intersection = linreg$coefficients[1]
      slope=linreg$coefficients[2]
      sum=summary(linreg)
      pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
      adjr2=sum$adj.r.squared
      if(m <= 1 || step >= m){
        adjr2s=c(adjr2s, adjr2)
        aics=c(aics,AIC(linreg))
        plottime=c(plottime,time[step])
      }
    }
    if(plot == TRUE){
      if(aic == TRUE){
        plot(plottime,aics,ylim=c(0,max(aics)),xlab="Time step",ylab="AIC",main=paste("Prediction accuracy",header), pch="+",type="b",col="blue")
      }else{
        plot(plottime,adjr2s,ylim=c(0,1),xlab="Time step",ylab="Adjusted R2",main=paste("Prediction accuracy",header), pch="+",type="b",col="blue")
      }
    }
  }
}
