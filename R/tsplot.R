#' @title Time Series Plot
#' @description Plot the time series row-wise.
#' @param x the matrix of time series
#' @param time.given if true, then the column names are supposed to hold the time units
#' @param num the number of rows to plot (starting from the first row)
#' @param sample.points sample.points
#' @param mode lines (default), pcoa (a PCoA plot with arrows showing the community trajectory) or bars (a stacked barplot for each sample)
#' @param dist the distance to use for the PCoA plot
#' @param identifyPoints click at points in the PCoA plot to identify them (using function identify), not active when noLabels is TRUE
#' @param topN number of top taxa to be plotted for mode bars
#' @param groups group membership vector with as many entries as samples, only supported for mode bars and pcoa
#' @param hideGroups compute PCoA with all data, but do not show members of selected groups; expects one integer per group and consistency with groups parameter, only supported for mode pcoa
#' @param legend add a legend
#' @param header string added to the plot title
#' @param labels use the provided labels in the PCoA plot
#' @param noLabels do not use any labels in the PCoA plot
#' @param perturb a perturbation object (adds polygons in mode lines highlighting the perturbation periods and colors dots in the PCoA plot)
#' @param \\dots Additional arguments passed to plot()
#' @examples
#' N=50
#' A=modifyA(generateA(N, c=0.1, d=-1),perc=70,strength="uniform",mode="negpercent")
#' out.ricker=ricker(N,A=A,y=generateAbundances(N,mode=5,prob=TRUE),K=rep(0.1,N), sigma=-1,tend=500)
#' tsplot(out.ricker)
#' tsplot(out.ricker[,1:20],mode="bars",legend=TRUE)
#' tsplot(out.ricker[,1:50],mode="pcoa")
#' @export
tsplot <- function(x, time.given=FALSE, num=nrow(x), sample.points=c(), mode="lines", dist="bray", identifyPoints=FALSE, topN=10, groups=c(), hideGroups=c(), legend=FALSE, header="", labels=c(), noLabels=FALSE, perturb=NULL, ...){
  col.vec = seq(0,1,1/nrow(x))
  my.colors = hsv(col.vec)
  main=paste("Community time series",header)
  xlab="Time points"
  ylab="Abundance"
  if(is.null(rownames(x))){
    rownames=c()
    for(row.index in 1:nrow(x)){
      rownames=c(rownames, row.index)
    }
    rownames(x)=rownames
  }
  if(time.given == TRUE){
    time=as.numeric(colnames(x))
  }else{
    time=c(1:ncol(x))
  }
  if(mode=="lines"){
    plot(time,as.numeric(x[1,]),ylim = range(x, na.rm = T),xlab = xlab, ylab = ylab, main = main, type = "l", col = my.colors[1], ...)
    # loop over rows in data
    for(i in 2:num){
      lines(time,as.numeric(x[i,]), col = my.colors[i], type = "l", ...)
    }
    # if non-empty, loop over sample points
    if(length(sample.points)>0){
      for(i in 1:length(sample.points)){
        abline(v=sample.points[i],col="gray", ...)
      }
    }
    if(!is.null(perturb)){
      for(perturbCounter in 1:length(perturb$times)){
        if(perturb$durations[perturbCounter]==1){
          abline(v=perturb$times[perturbCounter], col="red")
        }else{
          startTimePerturb=perturb$times[perturbCounter]
          stopTimePerturb=startTimePerturb+perturb$durations[perturbCounter]
          xx=c(rep(startTimePerturb,2),rep(stopTimePerturb,2))
          yy=c(0,rep(max(x,na.rm=TRUE),2),0)
          col1=rgb(1,0,0,0.3)
          polygon(xx,yy,col=col1)
        }
      }
    }
    if(legend == TRUE){
      legend("right",legend=rownames(x), lty = rep(1,nrow(x)), col = my.colors, merge = TRUE, bg = "white", text.col="black")
    }
  }else if(mode=="pcoa"){
    pcoa.res=vegan::capscale(data.frame(t(x))~1,distance=dist)
    colors=c()
    colors.copy=c()
    arrowColors=c()
    if(!is.null(perturb)){
      perturbCounter=1
      durationCounter=1
      perturbationOn=FALSE
      perturbIndices=c()
      defaultColor=rgb(0,1,0,0.5)
      perturbColor=rgb(1,0,0,0.5)
      for(timepoint in 1:ncol(x)){
        applied=applyPerturbation(perturb=perturb,t=timepoint, perturbCounter=perturbCounter, durationCounter=durationCounter, perturbationOn=perturbationOn, ori.growthrates = c(), abundances=c())
        durationCounter=applied$durationCounter
        perturbCounter=applied$perturbCounter
        perturbationOn=applied$perturbationOn
        if(perturbationOn==TRUE){
          colors=c(colors,perturbColor)
          perturbIndices=c(perturbIndices,timepoint)
        }else{
          colors=c(colors,defaultColor)
        }
      }
      # color preceding arrows also
      arrowColors=colors
      for(perturbIndex in perturbIndices){
        if(perturbIndex>0){
          arrowColors[(perturbIndex-1)]=perturbColor
        }
        if(perturbIndex<ncol(x)){
          nextIndex=perturbIndex+1
          if(length(which(perturbIndices==nextIndex))==0){
            arrowColors[perturbIndex]=defaultColor
          }
        }
      }
    }else{
      colors="black"
      arrowColors=rep("black",ncol(x))
    }
    # xlim with margin for the legend
    if(legend==TRUE){
      xlim=c(min(pcoa.res$CA$u[,1]),max(pcoa.res$CA$u[,1]+0.1))
    }else{
      xlim=c(min(pcoa.res$CA$u[,1]),max(pcoa.res$CA$u[,1]))
    }
    ylim=c(min(pcoa.res$CA$u[,2]),max(pcoa.res$CA$u[,2]))
    #labelNames=c(1:ncol(x))
    labelNames=time
    groups.copy=groups
    # color points according to groups
    if(length(groups)>0){
      groupNum=length(unique(groups))
      #print(paste("group number: ",groupNum))
      col.vec = seq(0,1,1/groupNum)
      hues = hsv(col.vec)
      hueCounter=1
      prevGroup=groups[1]
      colors=c()
      for(group.index in 2:length(groups)){
        colors=c(colors,hues[hueCounter])
        if(prevGroup!=groups[group.index]){
          hueCounter=hueCounter+1
        }
        prevGroup=groups[group.index]
      }
      colors.copy=colors
      #print(colors)
      if(length(hideGroups)>0){
        hidden.sample.indices=c()
        for(hidden.group.index in hideGroups){
          hidden.sample.indices=c(hidden.sample.indices, which(groups==hidden.group.index))
        }
        visible.sample.indices=setdiff(1:length(groups),hidden.sample.indices)
        pcoa.res$CA$u=pcoa.res$CA$u[visible.sample.indices,]
        labelNames=labelNames[visible.sample.indices]
        groups.copy=groups.copy[visible.sample.indices]
        colors.copy=colors[visible.sample.indices]
        colors.copy[length(colors.copy)]=colors.copy[length(colors.copy)-1]
        if(length(labels)> 0){
          labels=labels[visible.sample.indices]
        }
      }
    } # groups provided
    plot(pcoa.res$CA$u[,1:2], main=paste("PCoA",header), xlim=xlim, ylim=ylim, xlab="PCoA1", ylab="PCoA2", col=colors.copy, ...)
    #points(x=pcoa.res$CA$u[nrow(pcoa.res$CA$u),1],y=pcoa.res$CA$u[nrow(pcoa.res$CA$u),2], col=colors.copy[length(colors.copy)-1])
    for(i in 1:(nrow(pcoa.res$CA$u)-1)){
      if(length(groups)==0 || groups.copy[i]==groups.copy[i+1]){
        arrows(x0=pcoa.res$CA$u[i,1],y0=pcoa.res$CA$u[i,2],x1=pcoa.res$CA$u[i+1,1],y1=pcoa.res$CA$u[i+1,2], length=0.08, col=arrowColors[i])
      }
    }
    #print(pcoa.res$CA$u[nrow(pcoa.res$CA$u),1:2])
    if(noLabels==FALSE){
      if(length(labels)>0){
        labelNames=labels
        #print(labelNames)
      }
      text(pcoa.res$CA$u[,1],pcoa.res$CA$u[,2],labels=labelNames, pos=3, cex=0.8)
      if(identifyPoints==TRUE){
        identify(pcoa.res$CA$u[,1],pcoa.res$CA$u[,2])
      }
    }
    if(length(groups)>0 && legend==TRUE){
      legend("topright",legend=unique(groups),cex=0.8, pch = rep("-",length(unique(groups))), col = unique(colors), bg = "white", text.col="black")
    }
  }else if(mode=="bars"){

    if(is.null(colnames(x))){
      colnames=c()
      for(col.index in 1: ncol(x)){
        colnames=c(colnames,paste("T",col.index,sep=""))
      }
      colnames(x)=colnames
    }

    if(length(groups)>0){
      groupNum=unique(groups)
    }else{
      groupNum=1
    }

    # colors
    if(groupNum==1){
      colornumber=topN+1
    }else{
      # heuristic (different top taxa in different groups, maximal nrow colors needed)
      colornumber=nrow(x)/2
    }
    col.vec = seq(0,1,1/colornumber)
    my.colors = hsv(col.vec)
    color.index=1
    # gray color for non-top taxa
    colormap=list("Others"="#a9a9a9")

    # loop groups
    for(group.index in groupNum){
      if(length(groupNum) > 1){
        group.member.indices=which(groups==group.index)
        xsub=x[,group.member.indices]
        main=paste("Top",topN,"taxa of group",group.index,header)
      }else{
        xsub=x
        main=paste("Top",topN,"taxa",header)
      }
      rowsums=apply(xsub,1,sum)
      sorted=sort(rowsums,decreasing=TRUE,index.return=TRUE)
      sub.xsub=xsub[sorted$ix[1:topN],]
      misc=xsub[sorted$ix[(topN+1):nrow(xsub)],]
      if(ncol(xsub) > 1){
        misc.summed=apply(misc,2,sum)
        sub.xsub=rbind(sub.xsub,misc.summed)
      }else{
        misc.summed=sum(misc)
        sub.xsub=as.matrix(c(sub.xsub,misc.summed))
        rownames(sub.xsub)=c(rownames(xsub)[sorted$ix[1:topN]],"")
      }

      # add dummy columns for the legend (according to a heuristic)
      dummynum=round(ncol(xsub)/2)+round(ncol(xsub)/10)
      for(dummyindex in 1:dummynum){
        sub.xsub=cbind(sub.xsub,rep(NA,nrow(sub.xsub)))
        colnames(sub.xsub)[ncol(sub.xsub)]=""
      }
      rownames(sub.xsub)[nrow(sub.xsub)]="Others"
      colnames(sub.xsub)[(ncol(sub.xsub)-dummynum+1):ncol(sub.xsub)]=""

      # select colors
      selected.colors=c()
      for(name in rownames(sub.xsub)){
        if(name %in% names(colormap)){
          #print(paste("Found color",colormap[[name]],"for name",name))
          selected.colors=c(selected.colors,colormap[[name]])
        }else{
          selected.colors=c(selected.colors, my.colors[color.index])
          colormap[[name]]=my.colors[color.index]
          color.index=color.index+1
        }
      }
      barplot(as.matrix(sub.xsub),col=selected.colors, ylab="Abundance", main=main,cex.lab=0.8,las=2, ...)
      if(legend==TRUE){
        legend("topright",legend=rownames(sub.xsub),cex=0.9, bg = "white", text.col=selected.colors)
      }
    }
  }else{
    stop("Plot mode ",mode, "not supported. Supported modes are: lines, pcoa and bars")
  }
}
