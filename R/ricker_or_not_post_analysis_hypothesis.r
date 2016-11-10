######################################################
# COMPUTE COEFFICIENTs TO ASSESS THE NATURE OF THE TIMESERIES
######################################################
rm(list = ls())

library(ggplot2)
library(geigen)
library(MASS)

source("apply_limits.R")

#choose time series to analyze
timeseries.list= c(1:3,5:24,31:38,40:43) #1:3,5:24,31:43)
#choose sampling method (ie, choose max number of species kept, time skip, length kept)
Nmax= 60
tskip.list=c(1)       # j-loop  (if only c(1), we don't look for sub-sampling)
tend.max.list=c(3000) # m-loop  (if only one value, we don't look at the effect of shortening the ts)
j=1
m=1
wd<-getwd()
tableR=cbind(1,1,1,1,1)
list.hypo=c("a")
#loop over the different time series to analyze

for(ll in timeseries.list){
  print(ll)
  #load the desired time series 
  namefile=paste(ll,"_timeseries",sep="") #name of the folder where this matrix is stored 
  file.path= paste(wd,"/timeseries/",namefile,sep="") #full name of the path to this folder
  R<-t(read.table(paste(file.path,"/ts_",ll,"_complete.txt",sep="")))
  if(dim(R)[1]>500){
  R.no.transcient<-R[400:dim(R)[1],]
  meansR<-colMeans(R.no.transcient)
  temp<-t(matrix(rep(colMeans(R.no.transcient),dim(R.no.transcient)[1]),ncol=dim(R.no.transcient)[1]))
  Rnorm<-R.no.transcient-temp
  varianceR<-colMeans(Rnorm*Rnorm)
  varianceR<-sum(varianceR)/dim(R.no.transcient)[2]
  } else {
    varianceR<-NaN
    }
  namefile.specific = paste(ll,"_ts_N",Nmax,"_skip",tskip.list[j],'_tmax',tend.max.list[m],sep ="")
  file.path.specific = paste(file.path,"/",namefile.specific,sep ="")
  sd=dim(R)
  N=sd[2]
  tmax=sd[1]
  Rcoef<-read.table(paste(file.path.specific,"/",namefile,"_corr.txt",sep=''))
  Best<-read.table(paste(file.path.specific,"/",namefile,"_Best.txt",sep=''))
  
 
  load(paste(file.path.specific,"/mylistsize.Rda",sep=""))
  load(paste(file.path.specific,"/mylistbest1.Rda",sep=""))
  temp=list.numberkept.percOK[length(list.numberkept.percOK)]
  Z=temp[[1]][3]

  Nkept=Rcoef[1,]
  lastN.index=length(Nkept)
  slope=(Rcoef[7,lastN.index]-Rcoef[7,1])/(Rcoef[1,lastN.index]-Rcoef[1,1])*Nmax
 
  #write.matrix(schur.list,file=paste(namefile,"_schur.txt",sep=""))
  #write.matrix(R*Z,file=paste(namefile,"s",sep=''))
  #if(known.original.A==TRUE){
  #  write.matrix(K,file=paste(namefile,"_K.txt",sep=''))
  #  write.matrix(A,file=paste(namefile,"_A.txt",sep=''))
  #  write.matrix(corAB,file=paste(namefile,"_corAB.txt",sep=''))
  #}
  #write.matrix(Best,file=paste(namefile,"_Best.txt",sep=''))
  maxB=max(Best)/Z
  minB=min(Best)/Z
  deltaB=maxB-minB
  corrcoefm=Rcoef[7,lastN.index]
  

  
  
  
  if(varianceR< 1e-07){
    hypo<-'no dynamics'
  } else{
    if(deltaB>40){
      hypo<-'soc'
    } else if(deltaB<0.01){
      hypo<-'untb'
    } else if((deltaB>0.1) && (deltaB <6) && corrcoefm>.5){
      hypo<-'ricker'
    } else if((deltaB>0.01) && (deltaB <20) && corrcoefm<.1){
      hypo<-'dm'
    } else
      hypo<-'no hypo'
  }
  print(hypo)
  
  ## if deltaB is too big, you might to wrong guess, let look at that more closely:
  locmax=seq(1:length(list_BEST))
  locslope=seq(1:length(list_BEST))

  locprediction=rep('no hypo',length(list_BEST))
  if(deltaB>2000) {    # || hypo =='no hypo'
    for(ii in 1:length(list_BEST)){
      locmax[ii]<-max(list_BEST[ii][[1]]/Z)-min(list_BEST[ii][[1]]/Z)
      locslope[ii]<-(Rcoef[7,ii]-Rcoef[7,1])/(Rcoef[1,ii]-Rcoef[1,1])*Nmax
    }
    ## 
    for(jj in 1:length(list_BEST)){
    if(locmax[jj]>40 && locslope[jj]<0){
      locprediction[jj]<-'soc'
    }
    if(locmax[jj]<0.01){
      locprediction[jj]<-'untb'
    }
    if((locmax[jj]>0.1) && (locmax[jj] <6) && Rcoef[7,jj]>.5){
      locprediction[jj]<-'ricker'
    }
    if((locmax[jj]>0.01) && (locmax[jj]<20) && Rcoef[7,jj]<.1){
      locprediction[jj]<-'dm'
    }
    }
    #print(locprediction)
    ## 
  }
  
 
  ##
  
  
  
  
  
  if(locmax[length(list_BEST)]>2000){
    #print(locprediction)
    #print(length(list_BEST)-1)
    hypo<- locprediction[length(list_BEST)-1]
    deltaB<-max(list_BEST[length(list_BEST)-1][[1]]/Z)-min(list_BEST[length(list_BEST)-1][[1]]/Z)
      slope<- (Rcoef[7,lastN.index-1]-Rcoef[7,1])/(Rcoef[1,lastN.index-1]-Rcoef[1,1])*Nmax
      corrcoefm<-Rcoef[7,length(list_BEST)-1]
    print(hypo)
    if(locmax[length(list_BEST)-1]>1000){
      hypo->'no pred'
    }
  }
  
  list.hypo=append(list.hypo,hypo)
  tableR=rbind(tableR,cbind(ll,round(deltaB,digits=2),round(slope,digits=4),round(corrcoefm,digits=2),round(varianceR,digits=13)))
}

tableRdim1=dim(tableR)[1]
tableR=tableR[2:tableRdim1,]
colnames(tableR) <- c("ll","deltaB","slope","max corr","variance")
list.hypo<-list.hypo[2:length(list.hypo)]
rownames(tableR)<-list.hypo
