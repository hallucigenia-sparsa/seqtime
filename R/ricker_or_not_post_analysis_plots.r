################################################################################
# MAKE PLOTS OUT OF LIMITS ANALYSIS DONE WITH RICKER_OR_NOT_RICKER.R
################################################################################

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(cowplot)


transc=FALSE
full.length=FALSE

folder_name <-"/Users/sdebuyl/timeseries/"
name_file<-"_timeseries" #"timeseries"

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

################################################################################

# 1. Encode info about location of files with the results of the analysis
# 2. Loop over time series: goal is to put for each ts, all correlations coefs in a single data frame to do nice plots.

################################################################################


# 1. Encode info about location of files with the results of the analysis

#choose time series to analyze
timeseries.list=c(1) # 1,4,5,8,12:20,31:38,44:49,51:60
#choose sampling method (ie, choose max number of species kept, time skip, length kept)
Nmax= 60              # maximal number of species to include in the analysis 

subsambpling=FALSE   # set to TRUE to look only at time points [1501:1600]

mchosen= 1 # should be 1 (could be modified to visualized on plots the effect of changing the lenght of the time series) 
# tmaxV<-tend.max.list[mchosen]



# 2. Loop over time series

  for(ll in timeseries.list){
  
  index.plots.list<-1
  name_file_with_nbr=paste(ll,name_file,sep="") #name of the folder where this matrix is stored 
  file.path= paste(folder_name,name_file_with_nbr,sep="") #full name of the path to this folder 51_sliced_timeseries.txt
    
  R<-t(read.table(paste(file.path,'/',name_file_with_nbr,'.txt',sep="")))
  
  namefile="my_results"
  R[R<0]<-0
  sd=dim(R)
  N=sd[2]
  tmax=sd[1]
  tend=tmax
  tmaxV=tmax
  tskip.list=c(1)       # j-loop  (if we choose tskip.list=c(1), we don't look for sub-sampling)
  tend.max.list=c(tend) # m-loop  (if only one value, we don't look at the effect of shortening the ts)
  print(paste("N is ",N," and tmax is ",tmax,sep=""))
  
  # initiate lists of plots for the various time series. 
  index.plots.list<-1
  plots.list<-list()
  
  #initiate lists of autocorr et corr
  list.autocorr=list()
  list.autocorr2=list()
  list.autocorr3=list()
  list.autocorr4=list()
  list.corr=list()
  list.corr2=list()
  list.corr3=list()
  list.species.names=list()

  
  indexlist=1  # index for combining j & m loops. we will summarize results from different samplings in data.frames. 
  j=1
  m=1
  mchosen=1
  
  indexlist=1 
  namefile.specific = paste(ll,"_ts_N_",Nmax,"_skip_",tskip.list[j],'_tmax_',dim(R)[1],sep ="")
  file.path.specific = paste(file.path,"/",namefile.specific,sep ="")
    
  #extract correlation infos
        
  corr_table<-read.table(paste(file.path.specific,"/",namefile,"_corr.txt",sep=""))

  
  dfc1<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[2,])),sampling="correlation",tmax=toString(tend.max.list[m]),timestep=toString(1))
  list.corr[[indexlist]]<-dfc1
  dfc2<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[3,])),sampling="correlation",tmax=toString(tend.max.list[m]),timestep=toString(2))
  list.corr2[[indexlist]]<-dfc2
  dfc3<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[4,])),sampling="correlation",tmax=toString(tend.max.list[m]),timestep=toString(3))
  list.corr3[[indexlist]]<-dfc3
  
  df1<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[7,])),sampling="auto-correlation",tmax=toString(tend.max.list[m]),timestep=toString(1))
  list.autocorr[[indexlist]]<-df1
  df2<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[8,])),sampling="auto-correlation",tmax=toString(tend.max.list[m]),timestep=toString(2))
  list.autocorr2[[indexlist]]<-df2
  df3<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[9,])),sampling="auto-correlation",tmax=toString(tend.max.list[m]),timestep=toString(3))
  list.autocorr3[[indexlist]]<-df3
  
  #extract list of species kept      
        
  listspecieskeptname=paste(file.path.specific,"/namespecies.Rda",sep="")#list_BEST
  load(listspecieskeptname)
  
  list.species.names[[indexlist]]<-list.indices.OK[length(list.indices.OK)][[1]]
  
  #extract list of inferred matrices
  
  liste1name=paste(file.path.specific,"/mylistbest1.Rda",sep="")#list_BEST
  liste3name=paste(file.path.specific,"/mylistsize.Rda",sep="")#list.numberkept.percOK
  load(liste1name) 
  load(liste3name)
  list.numberkept.percOK=list.numberkept.percOK[!unlist(lapply(list_BEST, is.null))]
  list_BEST=list_BEST[!unlist(lapply(list_BEST, is.null))]
  
  
  ################################################################################  
  #### PLOT TABLE OF SPECIES KEPT
  ################################################################################  
  
  matrixspecies<-do.call(rbind,list.species.names)
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = .5)),
    colhead = list(fg_params=list(cex = 1.0)),
    rowhead = list(fg_params=list(cex = 1.0)))
  
  

  
  temp<-matrix(matrixspecies,ncol=10,byrow=TRUE)
  row_names=c("species 1-10")
  if(dim(matrixspecies)[2]>10){
    for(kkk in 2:(ceiling(dim(matrixspecies)[2]/10))){
      row_names<-c(row_names,paste("species ",kkk*10-9,"-",kkk*10,sep=""))
    }
  }
  rownames(temp)<-row_names
  
  pdf(paste(file.path.specific,"/",ll,"_most_abundant_species_Nmax_",Nmax,"_tmax_",tmaxV,".pdf",sep=""), width = 5, height =3)
  myt <- gridExtra::tableGrob(temp, theme = mytheme)
  plot.species.names<-grid.draw(myt)
  dev.off()
  
  ################################################################################
  #### PLOT BEST AS A FUNCTION OF N FOR A GIVEN SAMPLING #### 
  ################################################################################            
              
  listmatplots=list()
  indexmat=1

  listindexred<-seq(1,length(list_BEST),by=round(length(list_BEST)/3))

    for(k in listindexred){
       Z=list.numberkept.percOK[[k]][[3]]
       Best<-list_BEST[[k]]
       deltaB<-max(Best/Z)-min(Best/Z)
       scale.plot<-max(max(Best/Z),-min(Best/Z))
       listmatplots[[indexmat]]<-ggplot(melt(Best/Z), aes(Var1,Var2, fill=value)) + geom_tile()+ scale_fill_gradient2(low = "red", mid = "white", high = "blue",limits=c(-scale.plot, scale.plot)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_fixed() + labs(x = paste("Inferred interaction matrix for tmax = ",tend.max.list[mchosen]," \n and ", list.numberkept.percOK[[k]][[2]]," species kept (deltaB = ",signif(deltaB, digits = 2),")",sep=""), y= " ")
       indexmat<-indexmat+1
    }
  n <- length(listindexred)-1
  listmatplots[[indexmat]]<-grid.draw(myt)
  nCol <- 2
  pdf(paste(file.path.specific,"/",ll,"_inferred_interaction_matrix_Nmax_",Nmax,"_tmax_",tmaxV,".pdf",sep=""), width = 8, height = 12)
  do.call("grid.arrange", c(listmatplots,ncol=2))
  dev.off()
  
  
  indexlist<-indexlist+1

  
  ################################################################################
  #### PLOT CORRELATIONS #### 
  ################################################################################    
  
  plot.title.corr=paste(ll,"_corr_and_auto_corr_Nmax_",Nmax,"_tmax_",tmaxV,".png",sep="")
  
  totlist1auto=do.call(rbind.data.frame,list.autocorr)
  totlist1corr=do.call(rbind.data.frame,list.corr)
  
  totlist1=do.call(rbind.data.frame,list.autocorr)
  totlist2=do.call(rbind.data.frame,list.autocorr2)
  totlist3=do.call(rbind.data.frame,list.autocorr3)
  totlist4=do.call(rbind.data.frame,list.corr)
  totlist5=do.call(rbind.data.frame,list.corr2)
  totlist6=do.call(rbind.data.frame,list.corr3)
  totlist=do.call(rbind.data.frame,list(totlist1,totlist2,totlist3,totlist4,totlist5,totlist6))
  
  temp=subset(totlist,tmax==tmaxV)
  plots.list[[index.plots.list]]<-
  ggplot(temp, aes(x=x,y=y,colour=timestep, shape = sampling, group=interaction(timestep, sampling)))  + geom_line(size=.1) + geom_point()+ xlab("N") + ylab("(auto-)correlation coefficients")+ ggtitle(paste("Correlation coefficients for ",ll,"_timeseries",sep=""))+ theme(plot.title = element_text(size=11))+coord_cartesian(xlim = c(0, Nmax))+coord_cartesian(ylim = c(-1.1,1.1)) + scale_colour_discrete(name = "number of time points ahead",breaks=c("1","2","3"),labels=c("1","2","3")) +scale_shape_discrete(name  =" ")
  ggsave(paste(file.path.specific,"/",plot.title.corr,sep=""), width = 15, height = 15)
  
  
  dev.off()
  
  #### PLOT TOGETHER ALL INFO ABOUT TIME SERIES
  
  nCol <- 1
  NN<-min(min(dim(myt)[2],Nmax),30)
  pdf(paste(file.path.specific,"/",ll,"_summary_plot.pdf",sep=""), width = 8, height = 12)
  do.call("grid.arrange", c(c(listmatplots,grid.draw(myt[,1:NN]),plots.list[1]), ncol=nCol))
  dev.off()
  indexlist<-indexlist+1
  index.plots.list<-index.plots.list+2
} #end of ll loop


################################################################################
#### PLOT BEST AS A FUNCTION OF N FOR A GIVEN SAMPLING #### 
################################################################################            

listmatplots=list()
indexmat=1

listindexred<-seq(1,length(list_BEST),by=round(length(list_BEST)/3))

for(k in listindexred){
  Z=list.numberkept.percOK[[k]][[3]]
  Best<-list_BEST[[k]]
  deltaB<-max(Best/Z)-min(Best/Z)
  scale.plot<-max(max(Best/Z),-min(Best/Z))
  listmatplots[[indexmat]]<-ggplot(melt(Best/Z), aes(Var1,Var2, fill=value)) + geom_tile()+ scale_fill_gradient2(low = "red", mid = "white", high = "blue",limits=c(-scale.plot, scale.plot)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_fixed() + labs(x = paste("Inferred interaction matrix for tmax = ",tend.max.list[mchosen]," \n and ", list.numberkept.percOK[[k]][[2]]," species kept (deltaB = ",signif(deltaB, digits = 2),")",sep=""), y= " ")
  indexmat<-indexmat+1
}
n <- length(listindexred)-1
listmatplots[[indexmat]]<-grid.draw(myt)
nCol <- 2
pdf(paste(file.path.specific,"/",ll,"_inferred_interaction_matrix_Nmax_",Nmax,"_tmax_",tmaxV,".pdf",sep=""), width = 8, height = 12)
do.call("grid.arrange", c(listmatplots,ncol=2))
dev.off()


################################################################################
#### if known interaction matrix
################################################################################
#mat<-t(read.table(paste(file.path,'/',ll,'_interactionmatrix.txt',sep="")))
#sum(diag(cor(Ared,Best)))/sdred[2]




