##############################################################
# Apply Limits algorithm 
# 
# By : Sophie de Buyl, translated from a Mathematica code provided by from Charles Fisher and Pankaj Mehta. 
# Ref : Identifying Keystone Species in the Human Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression. PLoS One 9:e102451 (2014)
# Created: 29/01/2016
# 
# 
# input: needs to give an interaction matrix if you want to use Ricker model OR needs a time series
# output: infer the interaction matrix "Best" assuming an underlying Riker dynamics and compute correlation with original matrix if time series has been generated with a Ricker model 
# 
################################################################ 


apply_limits <- function(R, output.path,A=matrix(),known.original.A=FALSE,K=matrix(),percOKlist=c(.9,.8,.7,.6,.5,.4,.3,.2,.1,0),Nmax=120,tmax=5000,showplots=FALSE,remove.transcient=FALSE,allindices=seq(1,Nmax)){

wd<-getwd()
list_BEST=list()              # to store infered interaction matrices  (as a function of N)
list_A=list()                 # to store orignial interaction matrices (as a function of N)
list.numberkept.percOK=list()
list.indices.OK=list()         # keep list of indices corresponding to species kept. 

  
library(Matrix)
library(MASS)
library(matrixStats)
library(geigen)
library(lattice)
library(reshape2) 
library(ggplot2)
library(miscTools)
library(grid)
library(gridExtra)

######################################################
# manipulate data
######################################################

#remove absent species
indices.not.zero<-which(colMeans(R)!=0,arr.ind = T)
R<-R[,indices.not.zero]
sd<-dim(R)
tend<-sd[1]
N<-sd[2]


if(remove.transcient==TRUE){
if(tend>401){
  R<-R[100:tend,]
  sd<-dim(R)
  tend<-sd[1]
  N<-sd[2]
  select.time.points<- seq(1,tend-1,by=round(tend/300))
  R<-R[select.time.points,]
  #remove again absent species
  indices.not.zero=which(colMeans(R)!=0,arr.ind=T)
  R<-R[,indices.not.zero]
  sd<-dim(R)
  tend<-sd[1]
  N<-sd[2]
}
}

# mean abundancies 

meansR<-colMeans(R)
maxR<-max(meansR)

if(N>Nmax){
  sort.ind=sort(meansR, decreasing = TRUE, index.return = TRUE)$ix[1:Nmax]
  R<-R[,sort.ind]
  if(known.original.A==TRUE){
    A<-A[sort.ind,sort.ind]
  }
  sd<-dim(R)
  tend<-sd[1]
  N<-sd[2]
} else{
  sort.ind = indices.not.zero
}



meansR<-colMeans(R)
maxR<-max(meansR)




# normalise data
Z<-sum(R[tend,])
R<-R/Z


######################################################
# apply limits to decreasing number of species (removing the least abundant)
######################################################

# initialize loop
if(known.original.A==TRUE){
  corAB=0
}
cor.ts.estim.ts1<-0
cor.ts.estim.ts2<-0
cor.ts.estim.ts3<-0
cor.ts.estim.ts4<-0
cor.ts.estim.ts5<-0

autocorr1<-0
autocorr2<-0
autocorr3<-0
autocorr4<-0
autocorr5<-0

number.species.kept<-0
schur.list<-0
index=1
# loop over number of species kept
for(percOK in percOKlist){
  print(percOK)
  indicesOK<-which(meansR>percOK*maxR)
  
  if(length(indicesOK)>1 & is.element(length(indicesOK),number.species.kept)==FALSE){
    
    if(known.original.A==TRUE){
      Ared<-A[indicesOK,indicesOK]
    }
    
    Rred<-R[,indicesOK]
    sdred=dim(Rred)
    Kred<-colMeans(Rred)
    Zred<-sum(Rred[tend,])
    Rred<-Rred/Zred
    
    #compute estimation of the interaction matrix
    Best<-t(limits(Rred,1))
    for(i in 2:sdred[2]){
      Best<-cbind(Best,t(limits(Rred,i)))
    }
    Best<-t(Best)
  
    if(sum(is.na(Best))!=0){
      listind=which(is.na(Best),TRUE)
      species.to.remove=unique(as.vector(listind))
      indicesOK<-setdiff(indicesOK,species.to.remove)
      Best<-Best[indicesOK,indicesOK]
    }
    if(sum(is.na(Best))!=0){print(bug)}
                
                if(max(Re(eigen(Best)$values))>0){
                  max.Best<-max(Best)
                  Best<-remove.pos.eig(Best/max.Best)*max.Best#note: not sure why division by max.Best helps removing positive eigenvalues...
                  schur.list<-cbind(schur.list,percOK)
                }
                
                #compute interaction matrix old version. 
                sdred<-dim(as.matrix(Best))[1]
                # compute correlations between estimated time series and original one x points ahead 
                #
                # initiate loop
                y <- Rred[1,]
                RestimLoc1 <- Rred[1,]
                RestimLoc2 <- t(Rred[1:2,])
                RestimLoc3 <- t(Rred[1:3,])
                RestimLoc4 <- t(Rred[1:4,])
                RestimLoc5 <- t(Rred[1:5,])
              
                for(t in 2:tend){
                  tm1<- t-1
                 
                  RestimLoc1<-cbind(RestimLoc1,Rred[tm1,]*exp(Best%*%(matrix(Rred[tm1,]-Kred,ncol=1))))
                  RestimLoc2<-cbind(RestimLoc2,RestimLoc1[,t-1]*exp(Best%*%(RestimLoc1[,t-1]-Kred)))
                  RestimLoc3<-cbind(RestimLoc3,RestimLoc2[,t-1]*exp(Best%*%(RestimLoc2[,t-1]-Kred)))
                  RestimLoc4<-cbind(RestimLoc4,RestimLoc3[,t-1]*exp(Best%*%(RestimLoc3[,t-1]-Kred)))
                  RestimLoc5<-cbind(RestimLoc5,RestimLoc4[,t-1]*exp(Best%*%(RestimLoc4[,t-1]-Kred)))
                }
                # remove useless and transpose
                RestimLoc1<-t(RestimLoc1)
                RestimLoc2<-t(RestimLoc2[,1:tend])
                RestimLoc3<-t(RestimLoc3[,1:tend])
                RestimLoc4<-t(RestimLoc4[,1:tend])
                RestimLoc5<-t(RestimLoc5[,1:tend])
                # correlation shifted estimated time series & time series.
                cor.ts.estim.ts1<-cbind(cor.ts.estim.ts1,sum(diag(cor(Rred[2:tend,],RestimLoc1[2:tend,])))/length(indicesOK))
                cor.ts.estim.ts2<-cbind(cor.ts.estim.ts2,sum(diag(cor(Rred[3:tend,],RestimLoc2[3:tend,])))/length(indicesOK))
                cor.ts.estim.ts3<-cbind(cor.ts.estim.ts3,sum(diag(cor(Rred[4:tend,],RestimLoc3[4:tend,])))/length(indicesOK))
                cor.ts.estim.ts4<-cbind(cor.ts.estim.ts4,sum(diag(cor(Rred[5:tend,],RestimLoc4[5:tend,])))/length(indicesOK))
                cor.ts.estim.ts5<-cbind(cor.ts.estim.ts5,sum(diag(cor(Rred[6:tend,],RestimLoc5[6:tend,])))/length(indicesOK))
                # autocorrelation shifted time series.
                autocorr1<-cbind(autocorr1,sum(diag(cor(Rred[1:(tend-1),],Rred[2:tend,])))/length(indicesOK))
                autocorr2<-cbind(autocorr2,sum(diag(cor(Rred[1:(tend-2),],Rred[3:tend,])))/length(indicesOK))
                autocorr3<-cbind(autocorr3,sum(diag(cor(Rred[1:(tend-3),],Rred[4:tend,])))/length(indicesOK))
                autocorr4<-cbind(autocorr4,sum(diag(cor(Rred[1:(tend-4),],Rred[5:tend,])))/length(indicesOK))
                autocorr5<-cbind(autocorr5,sum(diag(cor(Rred[1:(tend-5),],Rred[6:tend,])))/length(indicesOK))
                
                if(known.original.A==TRUE){
                  corAB<-cbind(corAB,sum(diag(cor(Ared,Best)))/sdred[2])
                }
                
                
                
                
                ##########save Best for different number of species kept############
                
                #  if(showplots2==TRUE){
                #plot.title.AB=paste(output.path,"percOK",round(100*percOK),".png",sep="")
                
                if(known.original.A==TRUE){
                  A.temp=A[indicesOK,indicesOK]
                  list_BEST[[index]]=Best
                  list_A[[index]]=A.temp
                  list.numberkept.percOK[[index]]=c(percOK,length(indicesOK),Z)
                  list.indices.OK[[index]]=c(allindices[indices.not.zero[sort.ind[indicesOK]]])
                } else{
                  list_BEST[[index]]=Best
                  list.numberkept.percOK[[index]]=c(percOK,length(indicesOK),Z)
                  list.indices.OK[[index]]=c(allindices[indices.not.zero[sort.ind[indicesOK]]])
                }
                #}
                ##########end save Best for different number of species kept############
                number.species.kept<-cbind(number.species.kept,length(indicesOK))
                index<-index+1    
                
    #end of if condition
                
             } #end of if condition 
  
}#### end of loop over percentage. 



#to plot in "reanalyse.R" the infered interaction matrix as a function of number of species kept (and to know which species we kept.)
save(list_BEST,file=paste(output.path,"/mylistbest1.Rda",sep=""))
save(list_A,file=paste(output.path,"/mylistA2.Rda",sep=""))
save(list.numberkept.percOK,file=paste(output.path,"/mylistsize.Rda",sep=""))
save(list.indices.OK,file=paste(output.path,"/namespecies.Rda",sep="")) 
#remove first element of arrays (they have no meaning)

if(known.original.A==TRUE){
  corAB<-corAB[2:length(corAB)]
}
number.species.kept<-number.species.kept[2:length(number.species.kept)]
autocorr1<-autocorr1[-1]
autocorr2<-autocorr2[-1]
autocorr3<-autocorr3[-1]
autocorr4<-autocorr4[-1]
autocorr5<-autocorr5[-1]

cor.ts.estim.ts1<-cor.ts.estim.ts1[-1]
cor.ts.estim.ts2<-cor.ts.estim.ts2[-1]
cor.ts.estim.ts3<-cor.ts.estim.ts3[-1]
cor.ts.estim.ts4<-cor.ts.estim.ts4[-1]
cor.ts.estim.ts5<-cor.ts.estim.ts5[-1]


######################################################




######################################################
# save results in txt files
######################################################

setwd(output.path)

namefile="my_results"
write.matrix(schur.list,file=paste(namefile,"_schur.txt",sep=""))
write.matrix(R*Z,file=paste(namefile,"s",sep=''))
if(known.original.A==TRUE){
  write.matrix(K,file=paste(namefile,"_K.txt",sep=''))
  write.matrix(A,file=paste(namefile,"_A.txt",sep=''))
  write.matrix(corAB,file=paste(namefile,"_corAB.txt",sep=''))
}
write.matrix(Best,file=paste(namefile,"_Best.txt",sep=''))
write.matrix(rbind(number.species.kept,cor.ts.estim.ts1,cor.ts.estim.ts2,cor.ts.estim.ts3,cor.ts.estim.ts4,cor.ts.estim.ts5,autocorr1,autocorr2,autocorr3,autocorr4,autocorr5),file=paste(namefile,"_corr.txt",sep=""))
write.matrix(c(allindices[indices.not.zero[sort.ind[indicesOK]]]),file=paste(namefile,"_name_species_kept.txt",sep=''))

######################################################
# PLOTS
######################################################

plot.title.AB<-paste(namefile,"_AB.png",sep="")
plot.title.corrAB<-paste(namefile,"_corr_AB.png",sep="")
plot.title.ts<-paste(namefile,"_ts.png",sep="")
plot.title.AB<-paste(namefile,"_AB.png",sep="")
plot.title.corr<-paste(namefile,"_corr_coefs.png",sep="")
plot.title.ts.loc<-paste(namefile,"_loc_pred.png",sep="")


######################################################
# plot corAB
######################################################

if(showplots==TRUE){
if(known.original.A==TRUE){
  
  
  
        png(plot.title.corrAB)
        plot(number.species.kept,corAB,main="Correlation coefficients \n between inferred and original interaction matrices",xlab="number of species kept",ylab="correlation coefficient",xlim=c(0,N),ylim=c(-1,1),axes=FALSE)
        axis(1,at=1:N)
        axis(2,at=seq(-1,1,by=.2))
        grid()
        dev.off()
  
        
        

  scale.plot<-max(c(max(A),-min(A),max(Best/Z),-min(Best/Z)))
  
  p1<-ggplot(melt(A), aes(Var1,Var2, fill=value)) + geom_raster()+ scale_fill_gradient2(low = "red", mid = "white", high = "blue",limits=c(-scale.plot, scale.plot)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_fixed() + coord_fixed() + labs(x = "Original interaction matrix", y= " ")
  p2<-ggplot(melt(Best/Z), aes(Var1,Var2, fill=value),xlab="Inferred interaction matrix",ylab="y") + geom_raster()+ scale_fill_gradient2(low = "red", mid = "white", high = "blue",limits=c(-scale.plot, scale.plot))+  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_fixed() + coord_fixed()+ labs(x = "Inferred interaction matrix", y= " ") 
  ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
  ggsave(plot.title.AB,grid.arrange(p1,p2,ncol=2))
  
} else{
  scale.plot<-max(c(max(Best/Z),-min(Best/Z)))
  ggplot(melt(Best/Z), aes(Var1,Var2, fill=value),xlab="Inferred interaction matrix",ylab="y") + geom_raster()+ scale_fill_gradient2(low = "red", mid = "white", high = "blue",limits=c(-scale.plot, scale.plot))+  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_fixed() + coord_fixed()+ labs(x = "Inferred interaction matrix", y= " ") 
  ggsave(plot.title.AB)
  }


######################################################
# plots corr ts
######################################################
my.colors<-c("#FF0000FF","#0066FFFF","#CC00FFFF","#00FF66FF","#CCFF00FF")



test_corr <- data.frame(
  corrcoef1 = cor.ts.estim.ts1,
  corrcoef2 = cor.ts.estim.ts2,
  corrcoef3 = cor.ts.estim.ts3,
  corrcoef4 = cor.ts.estim.ts4,
  corrcoef5 = cor.ts.estim.ts5,
  autocor1 = autocorr1,
  autocor2 = autocorr2,
  autocor3 = autocorr3,
  autocor4 = autocorr4,
  autocor5 = autocorr5,
  date = number.species.kept#seq.Date(as.Date("2002-01-01"), by="1 month", length.out=100)
)


test_data_long <- melt(test_corr, id="date")  # convert to long format

#remove print
tempplot<-ggplot(data=test_data_long,ylim=c(-1,1),aes(x=date, y=value, colour=variable), xlab="Number species", ylab="corr")+
        geom_line()+geom_point(aes(color = variable,shape=variable))+xlim(c(0,N))+ylim(c(-1,1))+xlab("number of species")+ylab("correlation coefs")+ggtitle("Correlation coefficients") + scale_shape_manual(values=c(1,1,1,1,1,8,8,8,8,8))+scale_colour_manual(values=c(my.colors,my.colors))

ggsave(tempplot,file=plot.title.corr)


######################################################
# plot ts
######################################################




tsdf<- data.frame(
  speciesA = Rred[,1],
  date = seq(1,tend)
)

for(j in 2:N){
  tsdftemp<-data.frame(temp=Rred[,j],date=seq(1,tend))
  colnames(tsdftemp) <- c(paste("species",j,sep=""), "date")
  tsdf<-merge(tsdf,tsdftemp,by="date")
  
}

tsdf.melt <- melt(tsdf, id="date")  # convert to long format  ,by=c("ID","date")
#remove print
tempplot<-ggplot(data=tsdf.melt,ylim=c(0,max(Rred)),aes(x=date, y=value,colour=variable), xlab="Time", ylab="Abundancy")+
        geom_line(size=.1)+xlim(c(0,tend))+ylim(c(0,max(Rred)))+xlab("time")+ylab("abundancies")+ggtitle("Original Time Series")+theme(legend.position="none")

ggsave(tempplot,file=plot.title.ts)



######################################################
# plot loc predictions
######################################################




species.select <- sample(1:N)[1:3]
for(k in species.select){
  speciesN<-k
  name.loc.ts<-paste("loc.ts",toString(k),sep="")
  temp<-rep(NA,tend)
  for(m in seq(0,tend-5,by=round(tend/5)) ){
    start.time.point<-m
    temp[start.time.point:(start.time.point+5)]<-c(Rred[start.time.point,speciesN],RestimLoc1[start.time.point+1,speciesN],RestimLoc2[start.time.point+2,speciesN],RestimLoc3[start.time.point+3,speciesN],RestimLoc4[start.time.point+4,speciesN],RestimLoc5[start.time.point+5,speciesN])
  }
  assign(name.loc.ts,temp)
}





test_data <- data.frame(
  speciesA = Rred[,species.select[1]],
  speciesB = Rred[,species.select[2]],
  speciesC = Rred[,species.select[3]],
  loc_pred_A = get(paste("loc.ts",toString(species.select[1]),sep="")),
  loc_pred_B = get(paste("loc.ts",toString(species.select[2]),sep="")),
  loc_pred_C = get(paste("loc.ts",toString(species.select[3]),sep="")),
  date = seq(1,tend)#seq.Date(as.Date("2002-01-01"), by="1 month", length.out=100)
)


test_data_long <- melt(test_data, id="date")  # convert to long format

#remove print
tempplot<-ggplot(data=test_data_long,ylim=c(0,max(Rred)),aes(x=date, y=value, colour=variable), xlab="Time", ylab="Abundancy")+
  geom_line(size=.1)+geom_point(size=.5)+xlim(c(0,tend))+ylim(c(0,max(Rred)))+xlab("time")+ylab("abundancies")+ggtitle("3 species and local predictions")

ggsave(tempplot,file=plot.title.ts.loc)



}
setwd(wd)

######################################################
# write results in HTML format
######################################################
html.file<-paste(output.path,"/",namefile,"_html.txt",sep ='')
sink(html.file)  # open file to write data to - file will be deleted if it exists 

cat(paste("<tr> \n<td valign=", dQuote("top"),">",namefile,"</td>", "\n",collapse = ''))
cat("<td valign=", dQuote("top"),sprintf("> %d </td>\n",N))
cat("<td valign=",dQuote("top"),"><a href=",dQuote(sprintf("%s_ts.txt",namefile)),">R.txt</a> [<a href=",dQuote(plot.title.ts),"> png </a>] </td>\n")


if(known.original.A==TRUE){
  cat("<td valign=",dQuote("top"),"><a href=",dQuote(sprintf("%s_A.txt",namefile)),">A.txt</a> </td>\n")
  cat("<td valign=",dQuote("top"),"><a href=",dQuote(sprintf("%s_K.txt",namefile)),">K.txt</a> </td>\n")
} else{
  cat("<td valign=", dQuote("top"),"> not appl </td>\n")
  cat("<td valign=", dQuote("top"),"> not appl </td>\n") 
}
cat("<td valign=",dQuote("top"),"><a href=",dQuote(sprintf("%s_Best.txt",namefile)),"> Aestim.txt</a> [<a href=",dQuote(plot.title.AB),"> png </a>] </td>\n")
if(known.original.A==TRUE){
  cat("<td valign=",dQuote("top"),"><a href=",dQuote(plot.title.corrAB),"> plot.corrAB </a></td>\n")
  
} else{
  cat("<td valign=", dQuote("top"),"> not appl </td>\n")
  
}
cat("<td valign=",dQuote("top"),"><a href=",dQuote(plot.title.corr),"> plot.corr </a></td>\n")
cat("<td valign=",dQuote("top"),"><a href=",dQuote(plot.title.ts.loc),"> plot.ts.loc </a></td>\n")

cat("<td valign=",dQuote("top"),"><a href=",dQuote(sprintf("%s_corr.txt",namefile)),">corr.txt</a> </td>\n</tr> \n")
sink()                   # close file 

#should add a return function 
return(list(Best,rbind(number.species.kept,cor.ts.estim.ts1,cor.ts.estim.ts2,cor.ts.estim.ts3,cor.ts.estim.ts4,cor.ts.estim.ts5,autocorr1,autocorr2,autocorr3,autocorr4,autocorr5),c(allindices[indices.not.zero[sort.ind[indicesOK]]])))

}















################################################################ 
################################################################ 
################################################################ 
# Useful functions
################################################################ 
################################################################ 
################################################################ 






##############################################################
# Limits algorithm 
# 
# by Sophie de Buyl, translated from a Mathematica code provided by from Charles Fisher and Pankaj Mehta. 
# Ref : Identifying Keystone Species in the Human Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression. PLoS One 9:e102451 (2014)
# created: 29/01/2016
# 
# it provides an estimation of the "ith" line of a Ricker interaction matrix (apply code N time to get the full NxN matrix). 
# it requires the given of a time series
#
#
# Please adjust output.path before run!
#
# Should be use with apply_limits_script.R
#
# Needs packages:  matrixStats, MASS
################################################################ 

limits <- function(R,i){
  
  #### R is a time series with the columns corresponding to the species N, and 
  #### lines to different time points Ntp. If the first column of your time series
  #### i.e. dim(R)=Ntp N  
  #### labels time, remove it before applying limits (ie, do R<-R[:,2:end] ). 
  #### i is the column of the interaction matrix we try to reconstruct.
  
  #### output
  #### Beval : the output of the function is the estimation of the "ith" line of the interaction matrix 
  #### error : mean of the errors made on evalution of "y" using the estimated B matrix normalized by the variance of y ("one-time-step evalutaion"). 
  #### errorF : idem as error but without the normalization by the variance.
  listnumbkeysp<-c(); #list of number species kept. NOT MANDATORY
  
  # choices to be made:
  thresh <- .5; # orignially put to 5 (diminish to increase precision)
  r <- 100;     # the number of iterations used for Bagging 
  
  
  # manipulating R to put the data in the appropriated form to apply limits:
  
  R[R == 0] <-2.22507e-308 #we will have to take the log, so let's remove zero's. 2.220446e-16#
  sd<-dim(R)
  N<-sd[2]#number of species
  Ntp<-sd[1]-1#number of time points
  
  
  
  data<- R-t(kronecker(matrix(1,1,sd[1]),t(t(colMedians(R)))))#first N column of data matrix needed for limits
  data<-data[1:(sd[1]-1),]
  data<-cbind(data,(log(R[2:sd[1],i])-log(R[1:(sd[1]-1),i])))
  
  
  
  
  # variable initiation
  res<-array(0,c(N,1)) # array for storing results
  
  listspecies<-seq(1,N) # to construct the choices of excluded/includes species
  
  
  
  for(k in 1:r){ # we do r times the same thing to get r estimations the interaction matrix B
    
    #initialize covariates to be included/excluded from regression model
    if (i!= 1 & i!=N){
      c1<-1:(i-1)
      c2<-(i+1):N;
      excluded <-  c(c1,c2)
      included <- i;
    }
    
    if(i==1){
      excluded<-2:N
      included<-i
    }
    if(i==N){
      excluded<-1:(N-1)
      included<-N
    }
    
    #randomly partition data into training and test sets
    
    trainset <- as.matrix(sample(1:Ntp,round(Ntp/2)))
    testset <- as.matrix(setdiff(1:Ntp,trainset))
    
    data1 <- data[trainset,]
    data2 <- data[testset,]
    
    # loopN = 0
    # while(mean(as.matrix(data2[,dim(data2)[2]]))==0  && loopN != N-1)
    # {trainset <- as.matrix(sample(1:Ntp,round(Ntp/2)))
    # testset <- as.matrix(setdiff(1:Ntp,trainset))
    # 
    # data1 <- data[trainset,]
    # data2 <- data[testset,]
    # loopN<- loopN+1}
    
    test <- included
    results<- restrictedleastsquare(data1,data2,test) #perform the estimation
    errorEt <-results[[1]]
    B1t<-results[[2]]
    
    # loop that adds covariates(=species) as long as prediction error decreases by thresh
    xxx=1
    yyy=1
    if(N==2){
      B1temp=B1t
      }
    while(xxx == 1 && yyy != N-1){
      
      yyy=yyy+1
      
      #loop to create list of performances for all regression models including an additional covariate
      
      #initial the loop
      test <- c(included,excluded[1])
      errorlist<-restrictedleastsquare(data1,data2,test)[[1]]
      #the loop
      for(kk in 2:(length(excluded)-1)){
        test = c(included,excluded[kk])
        errorlist<-c(errorlist,restrictedleastsquare(data1,data2,test)[[1]])
      }
      
      
      #sort the list so that prev[[1]] has the lowest prediction error 
      ind<-which(errorlist == min(errorlist), arr.ind = TRUE)
      if(dim(as.matrix(ind))[1]>1){
        ind<-ind[[1]]
      }
      test <- c(included,excluded[ind])  #we choose the test set with the smallest error
      
      resulttemp <- restrictedleastsquare(data1, data2,test) #we re-run limits for that test set since we didn't store results
      errorEtemp <-resulttemp[[1]]
      B1temp<-resulttemp[[2]]
      # redefine included excluded covariates to account for including new covariate 
      included <- test
      excluded <- setdiff(union(listspecies,test),intersect(listspecies,test))
      
      # if improvement over the current best model by thresh, include new species, otherwise exit loop 
      
      condition.thresh<-100.0*(errorEtemp - errorEt)/errorEt> -thresh
      
      if( condition.thresh | is.nan(condition.thresh) | is.na(condition.thresh) ){ #we stop adding species
        xxx<-0
        listnumbkeysp<-c(listnumbkeysp,yyy)
      } else{ # we keep adding species
        errorEt <- errorEtemp# we update the error to be compared with. 
        B1t <- B1temp
        testt <- test 
      }
    } # end of while loop
    
    
    # store final regression coefficients in res 
    
    B <- matrix(0,N,1)
    B[test]<-t(B1temp)
    
    res<-cbind(res, B)
    
  } #end or r loop
  
  
  
  
  # Bagging step: output median of res 
  
  res<-res[,-1]
  if(dim(as.matrix(res))[2]>1){
    Beval<-t(as.matrix(rowMedians(res)))
  } else{
    Beval<-res # FMI dim(Beval)=1 6 (not sure I don't need a transposition here)
  }
  
  
  #listnumbkeysp 
  
  return(Beval)
}


##############################################################
# restrictedleastsquare
# function need for limits.r
# created: 29/01/2016
#
# by Sophie de Buyl, translated from a Mathematica code provided by from Charles Fisher and Pankaj Mehta. 
# Ref : Identifying Keystone Species in the Human Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression. PLoS One 9:e102451 (2014)
################################################################ 

restrictedleastsquare <- function(inbag,outbag,test){
  
  # This subroutine performs the linear regression and computes the error on the test set
  #"inbag" is a matrix containing the training set. One row of the matrix looks like {x_1(t) - u_1, ..., x_N(t) - 
  # u_N, ln x_i(t+1) - ln x_i(t) }, 
  #where u_i is the median abundance of species i. ;
  #"outbag" is a matrix containing the test set. 
  #It is in the same form as inbag.;
  #"test" is a vector of integers specifying which covariates are included in the regression. I.e. if species 1,2,3 are included then test = {1,2,3}.
  
  #perform linear regression on training set
  temp<-dim(inbag)
  lastcol<-temp[2]
  X1 <- as.matrix(inbag[,test])
  y1 <- as.matrix(inbag[,lastcol])
  B1 <- as.matrix(ginv(X1)%*%y1)#,tol=2e-308
  
  # calculate prediction error on test set 
  
  X2 <- as.matrix(outbag[,test])
  y2 <- as.matrix(outbag[,lastcol])
  errorE <- mean((y2-X2%*%B1)*(y2-X2%*%B1))/var(y2)
  result <- list(errorE,B1)
  return(result)
}

##############################################################
# remove.pos.eig 
# perform a schur decompostion of a matrix and flip positive real part of eigenvalues by adding a minus sign
# requires: geigen
# created: 29/01/2016
#
# by Sophie de Buyl
################################################################ 

remove.pos.eig<-function(A){

  sd<-dim(A)
  
  if(max(Re(eigen(A)$values))>0){
    diagt<-diag(sd[2])+0i
    
    schur.A<-gqz(A,diagt,"R")
    T<-schur.A$S%*%ginv(schur.A$T)
    rediag<-Re(diag(T))
    imdiag<-Im(diag(T))
    
    indicesP=rediag>0
    listind=1:sd[2]
    
    for(k in listind[indicesP]){
      T[k,k]<- complex(real=-Re(T[k,k]),imaginary=Im(T[k,k]))
    }
    
    A <- schur.A$Q %*% T %*% ginv(schur.A$Q)
    A<-Re(A)  
    return(A)
  }
}

  

#  adapted from Didier's code. 
#' Generate time series from Ricker model
#'
#' @param N number of species
#' @param A interaction matrix
#' @param K carrying capacity
#' @param y initial abundances
#' @param sigma noise level
#' @param tend number of time points
#' @param tskip number of initial time points to skip (to skip the transient)
#' @return a matrix with species abundances as rows and time points as columns
#' @examples
#' tsplot(ricker(10,generateA(10),K=rep(0.01,10)),type="l", header="ricker")
#' @seealso \code{\link{simuntb}} for the neutral model and \code{\link{glv}} for the generalized Lotka Volterra model
#' @export

ricker<-function(N, A, K=runif(N), y=runif(N), sigma=0.05, tend=300, tskip=0){
  bound=10^4 # above: explosion
  out=matrix(nrow=N, ncol=tend-tskip)
  out[,1]=y
  # simulate difference equation
  for(t in 2:tend){
    b=rlnorm(N,meanlog=0,sdlog=sigma)
    y=b*y*exp(A%*%(y-K))
    if(max(y) > bound){
      stop("Explosion")
    }
    else if(length(y[y<0]) > 0){
      stop("Species below 0!")
    }
    if(t > tskip){
      out[,t]=y
    }
  }
  
  return(list(out,K))
}

