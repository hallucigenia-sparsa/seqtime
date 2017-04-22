# Compare properties of time series. The function offers two comparisons: The first compares the initial taxon proportions with
# the proportions at the last time point (proportions are ensured by normalizing time series) and is carried out when the time
# series are provided. The second compares properties within one algorithm to test whether there is a significant dependency
# between two properties within time series generated with the same method. Note that for the first test, the recommmended method
# is diff, since correlations computed on even proportions will be NA.

# data: the result of compareTS
# timeseries (optional): can be collected with compareTS, when time series are given, the initial and final taxon proportions are compared
# algorithm: the algorithm for which the dependency between time series properties should be tested (ignored for timeseries)
# property1: the first of two properties for which the dependency should be tested (ignored for timeseries)
# property2: the second of two properties for which the dependency should be tested (ignored for timeseries)
# method: spearman, pearson, kendall or wilcoxon (wilcoxon only if at least one property is binary and not for timeseries), for timeseries also diff (sum of absolute differences between proportions)
#
# when timeseries are provided, the method returns the correlations or differences between initial and final taxon proportions for each experiment
#
# Examples:
# analyseTSComparison(table,algorithm="dm",property1 = "taylorslope",property2="initabundmode",method="wilcoxon")
# differences=analyseTSComparison(table, timeseries, method="diff")

analyseTSComparison<-function(data, timeseries=NULL, algorithm="ricker", property1="taylorr2", property2="sigma", method="spearman"){
  if(!is.data.frame(data)){
    data=as.data.frame(data)
  }
  if(!is.null(timeseries)){
    if(method=="wilcoxon"){
      stop("Wilcoxon not supported to compare time series.")
    }
    N=nrow(timeseries[["exp1"]])
    initProportionsBS=generateAbundances(N=N, mode=5, k=0.5, probabs=TRUE)
    initProportionsEven=rep(1/N,N)
    # initProportionsEven=initProportionsEven+runif(length(initProportionsEven),min=0,max=(min(initProportionsEven)/100))
    values=c()
    for(expId in 1:length(timeseries)){
      expName=paste("exp",expId,sep="")
      expTimeseries=normalize(timeseries[[expName]])
      lastTP=expTimeseries[,ncol(expTimeseries)]
      if(data$initabundmode[expId]==1){
        initProportions=initProportionsEven
      }else if(data$initabundmode[expId]==5){
        initProportions=initProportionsBS
      }else{
        stop(paste("Experiment",expId,"misses initial abundance mode!"))
      }
      if(method=="spearman" || method=="pearson" || method=="kendall"){
        value=cor.test(initProportions,lastTP,method=method)$estimate
      }else if(method=="diff"){
        value=sum(abs(initProportions-lastTP))
      }
      values=c(values,value)
    }
    return(values)
  }else{
    if(method=="diff"){
      stop("Method diff is only supported to compare time series.")
    }
    indices=which(data$algorithm==algorithm,arr.ind = TRUE)
    x=data[[property1]][indices]
    y=data[[property2]][indices]
    print(paste("Number of values:",length(indices)))
    if(method=="spearman" || method=="kendall" || method=="pearson"){
      cor.out=cor.test(x,y,method=method)
      print(paste("P-value:",cor.out$p.value))
      print(paste("Correlation:",cor.out$estimate))
    }else if(method=="wilcoxon"){
      valx=unique(x)
      valy=unique(y)
      if(length(valx)==2 && length(valy)==2){
        warning("Both variables can be treated as binary. One is selected randomly!")
      }
      if(length(valx)==2){
        categories=valx
        categoryAssignments=x
        values=y
      }else if(length(valy)==2){
        categories=valy
        categoryAssignments=y
        values=x
      }else{
        stop("Cannot apply Wilcoxon if there is not one variable with only two values!")
      }
      cat1=categories[1]
      cat2=categories[2]
      indicesCat1=which(categoryAssignments==cat1,arr.ind=TRUE)
      values1=values[indicesCat1]
      indicesCat2=which(categoryAssignments==cat2,arr.ind=TRUE)
      values2=values[indicesCat2]
      wilcox.out=wilcox.test(values1,values2)
      print(paste("P-value Wilcoxon",wilcox.out$p.value))
      print(paste("Category 1:",cat1))
      print("Values in category 1:")
      print(values1)
      print(paste("Category 2:",cat2))
      print("Values in category 2:")
      print(values2)
    }else{
      stop("Method not supported.")
    }
  }
}


