# Internal function: Collect fitting results to community time series.
#
# Different quality scores are computed for two fitting
# methods that infer interaction matrices, namely limits and soc.
# If the input settings are provided via the input.folder argument,
# the parameters obtained through fitting are compared to the
# real parameters.
# The following quality scores are computed:
# slope: The slope of the mean correlation of predicted and real time series versus species sub-set considered
# autoslope: The slope of the mean auto-correlation with lag 1 versus increasing species number, added in decreasing order of abundance
# corrAll: The mean cross-correlation between predicted and observed time series for all species considered
# maxcorr: The maximum among the mean cross-correlations between predicted and observed time series for different species sub-sets
# deltaAest: The range of values predicted for the estimated interaction matrx
# Acorr: The mean correlation of values in the predicted and the original interaction matrix (requires the input.folder).
# Note that each point of the predicted time series is obtained from the preceding point of the observed time series using Ricker with the inferred interaction matrix
#
# @param fit.folder the folder where fitting results are stored
# @param input.folder (optional) the folder where input settings for time series generated with function generateTS are stored
# param path.slices location of table that stores slice definitions (format: first column: start, second column: stop, NA: until end of time series, row: experiment identifier), if not provided, fit.method selected gives an error for soc
# param expIds the experiment identifiers of time series to be considered
# param dryRun skip time-consuming analyses: only deltaAest and Acorr are computed
# param setMissingIdsToNA if experimental identifiers have gaps, fill the gaps with NA values, e.g. c(1,3,4) will be filled to c(1,2,3,4) with missing values for experiment identifier 2
# param maxId to be used optionally in combination with setMissingIdsToNA (last identifier is not in the identifier set provided)
# param fit.method the fitting method (limits or soc)
# param fit.type the part of the time series used (all, selected, first), first only for limits
# param limits.round the round of LIMITs fitting (changes some folder names)
# return a table with goodness of fit scores, including slope, autoslope, corrAll, maxcorr, deltaAest and Acorr (Acorr is only computed if the input.folder is provided)

compareFit<-function(fit.folder="", input.folder="", path.slices="", expIds=c(), dryRun=FALSE, setMissingIdsToNA=FALSE, maxId=NA, fit.method="limits", fit.type="all", limits.round=2){
  if(fit.folder != ""){
    if(!file.exists(fit.folder)){
      stop(paste("The fitting result folder",fit.folder,"does not exist!"))
    }
    if(fit.method=="limits"){
      input.timeseries.folder=paste(fit.folder,"timeseries",sep="/")
      if(!file.exists(input.timeseries.folder)){
        stop("The fitting result folder does not have a time series subfolder!")
      }
    }else{
      input.timeseries.folder=fit.folder
    }
  }else{
    stop("Please provide the fitting result folder!")
  }

  slices=matrix()
  if(fit.method=="soc" && fit.type=="selected"){
    if(path.slices==""){
      stop("Please provide the slice definitions for SOC!")
    }else{
      slices=as.matrix(read.table(path.slices))
    }
  }

  if(input.folder != ""){
    if(!file.exists(input.folder)){
      stop(paste("The input folder",input.folder,"does not exist!"))
    }
    input.settings.folder=paste(input.folder,"settings",sep="/")
    if(!file.exists(input.settings.folder)){
      stop("The input folder does not have a settings subfolder!")
    }
  }

  limitsqualautocorslope=c()
  limitsqualslope=c()
  limitsqualmaxcorr=c()
  limitsqualcrosscorAll=c()
  limitsqualdeltaAest=c()
  limitsqualAcorr=c()

  gaps = c()
  loop = c()
  isGap = FALSE
  AestSOCName=""
  TSSOCName=""
  ExtinctionName=""
  ImmigrationName=""

  # identify gaps
  if(setMissingIdsToNA==TRUE){
    if(is.na(maxId)){
      maxExpId=max(expIds)
    }else{
      maxExpId=maxId
    }
    for(i in 1:maxExpId){
      loop=c(loop,i)
      if(i %in% expIds){
        # OK
      }else{
        gaps=c(gaps,i)
      }
    }
  }else{
    loop = expIds
  }

  for(expId in loop){

    print(paste("Processing identifier",expId))

    # is it a gap?
    if(expId %in% gaps){
      isGap=TRUE
    }else{
      isGap=FALSE
    }

    if(isGap==FALSE){

      if(fit.method=="limits"){
        input.timeseries.name=paste(expId,"timeseries",sep="_")

      }else{
        internalNumber=2
        if(expId>24){
          internalNumber=3
        }
        if(expId==42 || expId==43 || expId==44 || expId==46 || expId==47 || expId==49){
          internalNumber=2
        }
        if(fit.type == "selected"){
          input.timeseries.name=paste("SOI_TS_",expId,"_",internalNumber,"_subsample",sep="")
          AestSOCName=paste("soi_interaction_ts_",expId,"_",internalNumber,"_subsample.csv",sep="")
          TSSOCName=paste("soi_ts_",expId,"_",internalNumber,"_subsample.csv",sep="")
          ExtinctionName=paste("soi_extinction_ts_",expId,"_",internalNumber,"_subsample.csv",sep="")
          ImmigrationName=paste("soi_immigration_ts_",expId,"_",internalNumber,"_subsample.csv",sep="")
        }else if(fit.type=="all"){
          input.timeseries.name=paste("SOI_TS_",expId,sep="")
          AestSOCName=paste("soi_interaction_ts_",expId,".csv",sep="")
          TSSOCName=paste("soi_ts_",expId,".csv",sep="")
          ExtinctionName=paste("soi_extinction_ts_",expId,".csv",sep="")
          ImmigrationName=paste("soi_immigration_ts_",expId,".csv",sep="")
        }else{
          stop("Fit type first is not supported for soc.")
        }
      }
      input.timeseries.expId.folder=paste(input.timeseries.folder,input.timeseries.name,sep="/")
      print(input.timeseries.expId.folder)
      if(!file.exists(input.timeseries.expId.folder)){
        stop("The input time series folder does not have a subfolder for the input experiment identifier!")
      }

      # optionally read settings for current experiment
      if(input.folder != ""){
        input.settings.name=paste(expId,"settings",sep="_")
        input.settings.expId.folder=paste(input.settings.folder,input.settings.name,sep="/")
        if(!file.exists(input.settings.expId.folder)){
          stop("The input settings folder does not have a subfolder for the input experiment identifier!")
        }

        # read settings file
        input.settings.expId.file=paste(expId,"settings.txt",sep="_")
        settings.path=paste(input.settings.expId.folder,input.settings.expId.file,sep="/")
        if(!file.exists(settings.path)){
          stop(paste("The settings file",settings.path,"does not exist!"))
        }
        source(settings.path)

        # read input matrix
        if(Algorithm == "ricker" || Algorithm == "soc" || Algorithm == "glv"){
          source.expId = expId
          interactionmatrix.folder=input.settings.expId.folder
          # interaction matrix for current experiment was read from another experiment
          if(!is.na(Input_experiment_identifier)){
            source.expId=Input_experiment_identifier
            source.expId.folder=paste(source.expId,"settings",sep="_")
            interactionmatrix.folder=paste(input.settings.folder,source.expId.folder,sep="/")
          }
          A.name=paste(source.expId,"interactionmatrix.txt",sep="_")
          input.path.A=paste(interactionmatrix.folder,A.name,sep="/")
          print(paste("Reading interaction matrix from:",input.path.A,sep=" "))
          A=read.table(file=input.path.A,sep="\t",header=FALSE)
          A=as.matrix(A)
        }
      }

      # select fitting result
      if(fit.method=="limits"){
        subfolder = ""
        if((expId >= 25 && expId <= 30) || (expId >= 44 && expId <= 50)){
          subfolder=paste(expId,"ts_full_N60_skip1_tmax",sep="_")
          # special settings for second LIMITs round on slices
          if(limits.round==2){
            if(expId==25 || expId==28 || expId>=44){
              subfolder=paste(expId,"ts_N60_skip1_tmax100",sep="_")
            }else if(expId==26){
              subfolder=paste(expId,"ts_full_N60_skip1_tmax73",sep="_")
            }else if(expId==27){
              subfolder=paste(expId,"ts_full_N60_skip1_tmax37",sep="_")
            }else if(expId==29){
              subfolder=paste(expId,"ts_full_N60_skip1_tmax51",sep="_")
            }else if(expId==30){
              subfolder=paste(expId,"ts_full_N60_skip1_tmax26",sep="_")
            }
          }
        }else{
          if(fit.type=="all"){
            if(expId==56){
              subfolder=paste(expId,"ts_N60_skip1_tmax600",sep="_")
            }else{
              subfolder=paste(expId,"ts_N60_skip1_tmax3000",sep="_")
            }
          }else if(fit.type=="first"){
            # subfolder ID_ts_transcient_N60_skip_1_tmax100 is the first 100 time points
            subfolder=paste(expId,"ts_transcient_N60_skip1_tmax100",sep="_")
          }else if(fit.type=="selected"){
            subfolder=paste(expId,"ts_N60_skip1_tmax100",sep="_")
          }
        }
        expId.subfolder = paste(input.timeseries.expId.folder,subfolder,sep="/")
        expId.subfolder.corrfile=paste(expId.subfolder, paste(expId, "timeseries_corr.txt", sep="_"), sep="/")

        # format corrfile:
        # the first line labels the columns (number of species in ascending order)
        # the second the auto-correlation 1 step ahead
        # the 3rd, the auto-correlation 2 steps ahead
        # the 4th, the auto-correlation 3 steps ahead
        # the 5th, the auto-correlation 4 steps ahead
        # the 6th, the auto-correlation 5 steps ahead
        # the 7th, the cross-correlation 1 step ahead
        # the 8th, the cross-correlation 2 steps ahead
        # the 9th, the cross-correlation 3 steps ahead
        # the 10th, the cross-correlation 4 steps ahead
        # the 11th, the cross-correlation 5 steps ahead
        # Step ahead in cross-correlation: each time point of the predicted time series
        # was obtained from the preceding time point in the observed time series, but predicted
        # and observed time series are not shifted when computing their correlation
        corrs=read.table(file=expId.subfolder.corrfile,header=FALSE)
        corrs=as.matrix(corrs)
        specnum=as.numeric(corrs[1,])
        autocorr1=as.numeric(corrs[2,])
        #print(specnum)
        crosscorr1=as.numeric(corrs[7,])
        crosscorr1All=crosscorr1[length(crosscorr1)]
        expId.subfolder.Aestfile=paste(expId.subfolder, paste(expId, "timeseries_Best.txt", sep="_"), sep="/")
        Aest=read.table(file=expId.subfolder.Aestfile,header=FALSE)
        Aest=as.matrix(Aest)
        print(paste("Number of selected species",nrow(Aest)))
      }else{
        # read SOC fitting results
        expId.subfolder.Aestfile=paste(input.timeseries.expId.folder, AestSOCName, sep="/")
        print(expId.subfolder.Aestfile)
        Aest=read.table(file=expId.subfolder.Aestfile,header=TRUE)
        Aest=as.matrix(Aest)
        m.vector=read.table(file=paste(input.timeseries.expId.folder, ImmigrationName, sep="/"),skip=1,header=FALSE)
        m.vector=as.numeric(m.vector[,2])
        e.vector=read.table(file=paste(input.timeseries.expId.folder, ExtinctionName, sep="/"),skip=1,header=FALSE)
        e.vector=as.numeric(e.vector[,2])
        TS.location=paste(input.timeseries.expId.folder, TSSOCName, sep="/")
        oriTS=t(as.matrix(read.table(file=TS.location, header=TRUE)))
        if(fit.type=="selected"){
          start=slices[expId,1]
          end=slices[expId,2]
          if(is.na(end)){
            end=ncol(oriTS)
          }
          oriTS=oriTS[,start:end]
          print(paste("Time series sub-set length:",ncol(oriTS)))
        }
        if(dryRun==FALSE){
          limitsQualOut=limitsQuality(oriTS, Aest, norm=FALSE, type="soc", m.vector=m.vector, e.vector=e.vector, plot=FALSE, predict.stepwise=FALSE, noSchur=TRUE)
          specnum=limitsQualOut$taxonnum
          print(limitsQualOut$meancrosscor)
          crosscorr1=limitsQualOut$meancrosscor
          crosscorr1All=crosscorr1[length(crosscorr1)]
          autocorr1=limitsQualOut$meanautocor1
        }else{
          specnum=NA
          crosscorr1=NA
          crosscorr1All=NA
        }
      }

      # compute quality scores
      if(dryRun==FALSE){
        if(length(unique(crosscorr1))==1 && is.na(unique(crosscorr1))){
          linregslope=NA
        }else{
          #reg.data=data.frame(crosscorr1,specnum)
          linreg = lm(formula = crosscorr1~specnum)
          linregslope=linreg$coefficients[2]
        }
        if(length(unique(autocorr1))==1 && is.na(unique(autocorr1))){
          linregautoslope=NA
        }else{
          #reg.data.auto=data.frame(autocorr1,specnum)
          linreg.auto = lm(formula = autocorr1~specnum)
          linregautoslope=linreg.auto$coefficients[2]
        }
      }else{
        linregslope=NA
        linregautoslope=NA
      }
      deltaAest=max(Aest)-min(Aest)

      if(input.folder != ""){
        if(Algorithm == "glv" || Algorithm == "soc" || Algorithm == "ricker"){
          if(fit.method=="limits"){
            # get selected sub-set of A
            expId.subfolder.selectedSpecFile=paste(expId.subfolder,paste(expId,"timeseries_name_species_kept.txt", sep="_"),sep="/")
            selectedSpec=as.matrix(read.table(file=expId.subfolder.selectedSpecFile,header=FALSE))
            #print(selectedSpec[1:2,])
            A=A[selectedSpec[,1],selectedSpec[,1]]
          }
          Acorr=sum(diag(cor(A,Aest)))/nrow(Aest)
        }else{
          Acorr=NA
        }
        limitsqualAcorr=c(limitsqualAcorr,Acorr)
      }
      limitsqualmaxcorr=c(limitsqualmaxcorr,max(crosscorr1,na.rm=TRUE))
      limitsqualcrosscorAll=c(limitsqualcrosscorAll,crosscorr1All)
      limitsqualslope=c(limitsqualslope,linregslope)
      limitsqualdeltaAest=c(limitsqualdeltaAest,deltaAest)
      limitsqualautocorslope=c(limitsqualautocorslope, linregautoslope)
    }else{
      # treat gaps
      limitsqualcrosscorAll=c(limitsqualcrosscorAll,crosscorr1All)
      limitsqualmaxcorr=c(limitsqualmaxcorr,NA)
      limitsqualslope=c(limitsqualslope,NA)
      limitsqualdeltaAest=c(limitsqualdeltaAest,NA)
      limitsqualAcorr=c(limitsqualAcorr,NA)
      limitsqualautocorslope=c(limitsqualautocorslope,NA)
    }
  }

  # assemble table
  resulttable=list(limitsqualslope, limitsqualautocorslope,limitsqualcrosscorAll, limitsqualmaxcorr, limitsqualdeltaAest, limitsqualAcorr)
  names(resulttable)=c("slope","autoslope","corrAll","maxcorr","deltaAest", "Acorr")
  return(resulttable)
}
