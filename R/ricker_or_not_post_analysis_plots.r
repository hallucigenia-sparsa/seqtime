################################################################################
# MAKE PLOTS OUT OF LIMITS ANALYSIS DONE WITH RICKER_OR_NOT_RICKER.R
################################################################################

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

wd<-getwd()

################################################################################

# 1. Encode info about location of files with the results of the analysis
# 2. Loop over time series: goal is to put for each ts, all correlations coefs in a single data frame for do nice plots.

################################################################################




# 1. Encode info about location of files with the results of the analysis

timeseries.list=c(1:3,5:24,31:38,40:43) #
Nmax= 60 #the maximal number of species kept for the analysis
tskip.list=c(1)       # should be = c(1). we don't look for sub-sampling, to do that you need to add an extra loop (indexed by j in previous versions).
tend.max.list=c(3000) # should be only one number. if we want to look at different lenghts of the same time series, you need to add a loop (indexed by m in previous versions). 
mchosen= 1 # should be 1 (could be modified to visualized on plots the effect of changing the lenght of the time series) 
tmaxV<-tend.max.list[mchosen]



# 2. Loop over time series

  for(ll in timeseries.list){
  
  
  # initiate lists of plots for the various time series. 
  index.plots.list<-1
  plots.list<-list()
    
               
  
  #initiate lists of autocorr et corr
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
  
  # path to the ll time series
  namefile=paste(ll,"_timeseries",sep="") #name of the folder where this matrix is stored 
  file.path= paste(wd,"/timeseries/",namefile,sep="") #full name of the path to this folder
  
  # specify file for the given sampling method
  
  namefile.specific=paste(ll,"_ts_N",Nmax,"_skip",tskip.list[j],'_tmax',tend.max.list[m],sep="")
  file.path.specific=paste(file.path,"/",namefile.specific,sep="")
  file.path.ts=paste(file.path,"/ts_",ll,"_complete.txt",sep="")
  R<-read.table(file.path.ts)#ATTENTION A LA TRANsPOSITION   R<-read.table("/Users/sdebuyl/Dropbox/testhydracopy/timeseries/1_timeseries/1_timeseries_ts.txt")
  #R<-t(R)#Rtest<-read.table(paste(file.path,"/ts_1_complete.txt",sep=""))
  sd=dim(R)
  N=sd[2]
  tmax=sd[1]
  

        
  #extract correlation infos
        
  corr_table<-read.table(paste(file.path.specific,"/",namefile,"_corr.txt",sep=""))

  
  dfc1<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[7,])),sampling="correlation",tmax=toString(tend.max.list[m]),timestep=toString(1))
  list.corr[[indexlist]]<-dfc1
  dfc2<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[8,])),sampling="correlation",tmax=toString(tend.max.list[m]),timestep=toString(2))
  list.corr2[[indexlist]]<-dfc2
  dfc3<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[9,])),sampling="correlation",tmax=toString(tend.max.list[m]),timestep=toString(3))
  list.corr3[[indexlist]]<-dfc3
  
  df1<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[2,])),sampling="auto-correlation",tmax=toString(tend.max.list[m]),timestep=toString(1))
  list.autocorr[[indexlist]]<-df1
  df2<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[3,])),sampling="auto-correlation",tmax=toString(tend.max.list[m]),timestep=toString(2))
  list.autocorr2[[indexlist]]<-df2
  df3<-data.frame(x=as.numeric(unlist(corr_table[1,])),y=as.numeric(unlist(corr_table[4,])),sampling="auto-correlation",tmax=toString(tend.max.list[m]),timestep=toString(3))
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
  temp<-matrix(matrixspecies,ncol=10,byrow=TRUE)
  rownames(temp) <- c("species 1-10", "species 11-20","species 21-30","species 31-40","species 41-50","species 51-60")
  
  
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = .5)),
    colhead = list(fg_params=list(cex = 1.0)),
    rowhead = list(fg_params=list(cex = 1.0)))
  
  myt <- gridExtra::tableGrob(temp, theme = mytheme)
  pdf(paste(file.path.specific,"/",ll,"_most_abundant_species_Nmax_",Nmax,"_tend_",tmaxV,".pdf",sep=""), width = 5, height =3)
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
  pdf(paste(file.path.specific,"/",ll,"_inferred_interaction_matrix_Nmax_",Nmax,"_tend_",tmaxV,".pdf",sep=""), width = 8, height = 12)
  do.call("grid.arrange", c(listmatplots,ncol=2))
  dev.off()
  
  indexlist<-indexlist+1

  
  ################################################################################
  #### PLOT CORRELATIONS #### 
  ################################################################################    
  
  plot.title.corr=paste(ll,"_corr_and_auto_corr_Nmax_",Nmax,"_tend_",tmaxV,".png",sep="")
  
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
  plots.list[[index.plots.list+1]]<-
  ggplot(temp, aes(x=x,y=y,colour=timestep, shape = sampling, group=interaction(timestep, sampling)))  + geom_line(size=.1) + geom_point()+ xlab("N") + ylab("(auto-)correlation coefficients")+ ggtitle(paste("Correlation coefficients for ",ll,"_timeseries",sep=""))+ theme(plot.title = element_text(size=11))+coord_cartesian(xlim = c(0, Nmax))+coord_cartesian(ylim = c(-1.1,1.1)) + scale_colour_discrete(name = "number of time points ahead",breaks=c("1","2","3"),labels=c("1","2","3")) +scale_shape_discrete(name  =" ")
  ggsave(paste(file.path.specific,"/",plot.title.corr,sep=""), width = 5, height = 5)
 

} #end of ll loop









