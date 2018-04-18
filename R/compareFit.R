# Internal function: Collect fitting results to community time series.
#
# Different quality scores are computed for two fitting
# methods that infer interaction matrices, namely LIMITS and SOI fitting.
# If the input settings are provided via the input.folder argument,
# the parameters obtained through fitting are compared to the
# real parameters.
# The following quality scores are computed:
# slope: The slope of the mean correlation of predicted and real time series versus species sub-set considered
# autoslope: The slope of the mean auto-correlation with lag 1 versus increasing species number, added in decreasing order of abundance
# corrAll: The mean cross-correlation between predicted and observed time series for all species considered
# maxcorr: The maximum among the mean cross-correlations between predicted and observed time series for different species sub-sets
# deltaAest: The range of values predicted for the estimated interaction matrix
# Acorr: The mean correlation of values in the predicted and the original interaction matrix (requires the input.folder)
# Acon: The connectance of the inferred interaction matrix
# Apep: The positive edge percentage of the inferred interaction matrix
# Adim: The dimension of the inferred interaction matrix (i.e. the species number)
# Asum: The sum of absolute values in the inferred interaction matrix
# Arand: The mean of 10 iterations of mean correlation of known and random matrix (only reported if rand is set to true)
# Arandacc: Same as Arand, but instead of the correlation, the accuracy (mean of PPV and sensitivity) is reported (ignoring sign and strength; only reported if rand is set to true)
# oriAcon: The connectance of the original interaction matrix
# oriApep: The positive edge percentage of the original interaction matrix
# sens: sentivitity (the percentage of non-zero entries in the original matrix that are also in the inferred matrix), where sign and strength are ignored
# ppv: positive predictive value (the percentage of non-zero entries in the inferred matrix that are in the original matrix), where sign and strength are ignored
# acc: The accuracy of the inferred matrix (mean of PPV and sensitivity), where sign and strength are ignored
# Note that each point of the predicted time series is obtained from the preceding point of the observed time series using Ricker with the inferred interaction matrix, unless predict.stepwise is set to FALSE.
# The original time series are normalized before correlation is computed, but predicted time series are not normalized.
# The reason is that only a sub-set of the normalized time series is considered for prediction.
# Note also that random matrices do not preserve the diagonal, since LIMITS does not pre-set it. Inferred diagonal values are considered to be true positives.
#
# param fit.folder the folder where fitting results are stored
# param input.folder (optional) the folder where input settings and time series generated with function generateTS are stored
# param path.slices location of table that stores slice definitions (format: first column: start, second column: stop, NA: until end of time series, row: experiment identifier), if not provided, fit.method selected gives an error for soc
# param expIds the experiment identifiers of time series to be considered
# param rand compute the correlation of the known matrix with 10 randomly generated matrices of the same connectance as the known matrix and report the mean random correlation (only possible if input.folder is provided)
# param norm normalize time series (has only an impact if recomputeQual is TRUE)
# param sim similarity to use to assess discrepancy between observed and predicted time series (either r or kld, has only an impact if recomputeQual is TRUE)
# param recomputeQual if FALSE, quality scores involving time series prediction are read in from peviously generated files (limits) or are omitted (soc) - note that for limits, if no correlation file is found, the correlations are computed even if recomputeQual is false
# param predict.stepwise do a step-wise prediction (i.e. re-run model with result of previous step) instead of one run of the model (only relevant if recomputeQual is TRUE)
# param setMissingIdsToNA if experimental identifiers have gaps, fill the gaps with NA values, e.g. c(1,3,4) will be filled to c(1,2,3,4) with missing values for experiment identifier 2
# param maxId to be used optionally in combination with setMissingIdsToNA (last identifier is not in the identifier set provided)
# param fit.method the fitting method (limits or soc)
# param fit.type the part of the time series used (all, selected, first), first only for limits
# param limits.round the round of LIMITs fitting (changes some folder names)
# return a table with goodness of fit scores, including slope, autoslope, corrAll, maxcorr, deltaAest, Acorr, Acon, oriAcon, Apep, oriApep, Adim, Asum, Arand, percent (Acorr, oriAcon, oriApep & percent are only computed if the input.folder is provided, Arand in addition needs rand to be enabled)

compareFit<-function(fit.folder="", input.folder="", path.slices="", expIds=c(), rand=FALSE, norm=TRUE, sim="r", recomputeQual=FALSE, predict.stepwise=TRUE, setMissingIdsToNA=FALSE, maxId=NA, fit.method="limits", fit.type="all", limits.round=2){
  if(fit.folder != ""){
    if(!file.exists(fit.folder)){
      stop(paste("The fitting result folder",fit.folder,"does not exist!"))
    }
    if(fit.method=="limits"){
      input.timeseries.folder=paste(fit.folder,"timeseries",sep="/")
      if(fit.type=="noise"){
        input.timeseries.folder=paste(fit.folder,"Poisson",sep="/")
      }
      if(!file.exists(input.timeseries.folder)){
        stop("The fitting result folder does not have a time series subfolder!")
      }
    }else{
      input.timeseries.folder=fit.folder
    }
  }else{
    stop("Please provide the fitting result folder!")
  }

  if(fit.method=="soc" && predict.stepwise==TRUE){
    warning("Step-wise prediction not supported for SOI. It is set to false.")
    predict.stepwise=FALSE
  }

  if(recomputeQual == TRUE && input.folder==""){
    stop("To recompute the quality, the input folder is needed!")
  }

  if(fit.method == "soc" && recomputeQual == FALSE && input.folder==""){
    stop("To compute SOC fitting quality scores without recomputing time series, the input folder is needed for matrix quality computation!")
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
    input.ori.timeseries.folder=paste(input.folder,"timeseries",sep="/")
    print(paste("Observed time series folder:",input.ori.timeseries.folder))
    if(!file.exists(input.ori.timeseries.folder)){
      stop("The input folder does not have a timeseries subfolder with original time series!")
    }
  }

  randiter=10

  limitsqualautocorslope=c()
  limitsqualslope=c()
  limitsqualmaxcorr=c()
  limitsqualcrosscorAll=c()
  limitsqualdeltaAest=c()
  limitsqualAcorr=c()
  limitsqualConn=c()
  limitsqualPep=c()
  limitsqualdimA=c()
  limitsqualsumA=c()
  limitsqualsrandA=c()
  limitsqualsrandAperc=c()
  limitsqualAcc=c()
  limitsqualPPV=c()
  limitsqualSens=c()

  Acorr=NA
  oriACon=NA
  oriAPep=NA
  Arand=NA
  Arandperc=NA
  acc=NA
  sens=NA
  ppv=NA

  oriConn=c()
  oriPep=c()

  gaps = c()
  loop = c()
  isGap = FALSE
  AestSOCName=""
  TSSOCName=""
  ExtinctionName=""
  ImmigrationName=""

  Algorithm=""
  Input_experiment_identifier=NA
  Sampling_frequency=NA
  init_abundance_mode=NA
  theta=NA
  immigration_rate_Hubbell=NA
  deathrate_Hubbell=NA


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

  print(gaps)

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
        if(fit.type=="noise"){
          input.timeseries.name=paste(expId,"pois","timeseries",sep="_")
        }
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

        input.ori.timeseries.name=paste(expId,"timeseries",sep="_")
        input.ori.timeseries.expId.folder=paste(input.ori.timeseries.folder,input.ori.timeseries.name,sep="/")
        if(!file.exists(input.ori.timeseries.expId.folder)){
          stop("The input timeseries folder with original time series does not have a subfolder for the input experiment identifier!")
        }

        # read settings file
        input.settings.expId.file=paste(expId,"settings.txt",sep="_")
        settings.path=paste(input.settings.expId.folder,input.settings.expId.file,sep="/")
        if(!file.exists(settings.path)){
          stop(paste("The settings file",settings.path,"does not exist!"))
        }
        print(paste("Reading settings file from:",settings.path))
        source(settings.path, local=TRUE)
        #print(paste("Algorithm",Algorithm))
        #print(paste("steps",steps))

        # read time series
        input.ori.timeseries.expId.file=paste(expId,"timeseries.txt",sep="_")
        ori.timeseries.path=paste(input.ori.timeseries.expId.folder,input.ori.timeseries.expId.file,sep="/")
        if(!file.exists(ori.timeseries.path)){
          stop(paste("The timeseries file",ori.timeseries.path,"does not exist!"))
        }

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
          print(paste("Reading known interaction matrix from:",input.path.A,sep=" "))
          A=read.table(file=input.path.A,sep="\t",header=FALSE)
          A=as.matrix(A)
        }
      }

      # select fitting result
      if(fit.method=="limits"){
        subfolder = ""
        if(fit.type=="noise"){
          subfolder=paste(expId,"ts_N_60_skip_1_tmax_3000",sep="_")
          if(expId==45 || expId==46 || expId==47 || expId==48){
            subfolder=paste(expId,"ts_N_60_skip_1_tmax_300",sep="_")
          }
          if(expId==57){
            subfolder=paste(expId,"ts_N_60_skip_1_tmax_600",sep="_")
          }
        }else if(fit.type=="last"){
          # nothing
        }else{
          if((expId >= 25 && expId <= 30) || (expId >= 44 && expId <= 50 && fit.type!="first")){
            subfolder=paste(expId,"ts_full_N60_skip1_tmax",sep="_")
            # special settings for second LIMITs round on slices
            if(limits.round==2){
              if(expId==25 || expId==28){
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
        }
        expId.subfolder = paste(input.timeseries.expId.folder,subfolder,sep="/")
        expId.subfolder.corrfile=paste(expId.subfolder, paste(expId, "timeseries_corr.txt", sep="_"), sep="/")
        qualType="ricker"
        m.vector=c()
        e.vector=c()
        expId.subfolder.selectedSpecFile=paste(expId.subfolder,paste(expId,"timeseries_name_species_kept.txt", sep="_"),sep="/")
        if(fit.type=="noise" || fit.type=="last"){
          expId.subfolder.selectedSpecFile=paste(expId.subfolder,paste("my_results_name_species_kept.txt"),sep="/")
          expId.subfolder.corrfile=paste(expId.subfolder, paste("my_results_corr.txt"), sep="/")
        }
        print(expId.subfolder.selectedSpecFile)
        selectedSpec=as.matrix(read.table(file=expId.subfolder.selectedSpecFile,header=FALSE))
        if(recomputeQual==FALSE && file.exists(expId.subfolder.corrfile)){
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
        }else{
          oriTS=as.matrix(read.table(file=ori.timeseries.path))
          if(norm==TRUE){
            # normalize before discarding taxa
            oriTS=normalize(oriTS)
          }
          # only keep selected species
          oriTS=oriTS[selectedSpec[,1],]
          print(paste("Reading original time series from file:",ori.timeseries.path))
        }
        expId.subfolder.Aestfile=paste(expId.subfolder, paste(expId, "timeseries_Best.txt", sep="_"), sep="/")
        if(fit.type=="noise" || fit.type=="last"){
          expId.subfolder.Aestfile=paste(expId.subfolder, paste("my_results_Best.txt"), sep="/")
        }
        print(paste("Reading predicted interaction matrix from:",expId.subfolder.Aestfile))
        Aest=read.table(file=expId.subfolder.Aestfile,header=FALSE)
        Aest=as.matrix(Aest)
        AestDim=ncol(Aest)
        Acon=getConnectance(Aest)
        Apep=getPep(Aest)
        Asum=sum(abs(Aest))
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
        if(norm==TRUE){
          oriTS=normalize(oriTS)
        }
        qualType="soc"
        if(recomputeQual==FALSE){
          specnum=NA
          crosscorr1=NA
          crosscorr1All=NA
        }
      }

      if(recomputeQual==TRUE && input.folder != ""){

        # select sub-set
        if(fit.type=="selected"){
          start=slices[expId,1]
          end=slices[expId,2]
          if(is.na(end)){
            end=ncol(oriTS)
          }
          oriTS=oriTS[,start:end]
          print(paste("Time series sub-set length:",ncol(oriTS)))
        }

        # normalizing is done before, because for LIMITS, some taxa are discarded
        limitsQualOut=limitsQuality(oriTS, Aest, norm=FALSE, sim=sim, type=qualType, m.vector=m.vector, e.vector=e.vector, plot=FALSE, predict.stepwise=predict.stepwise, noSchur=TRUE)
        specnum=limitsQualOut$taxonnum
        print("Agreement observed and predicted time series:")
        print(limitsQualOut$meancrosscor)
        crosscorr1=limitsQualOut$meancrosscor
        crosscorr1All=crosscorr1[length(crosscorr1)]
        autocorr1=limitsQualOut$meanautocor1
      }

      # compute quality scores
      if(fit.method=="limits" || (fit.method=="soc" && recomputeQual==TRUE)){
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
        print(paste("Algorithm",Algorithm))
        if(Algorithm == "glv" || Algorithm == "soc" || Algorithm == "ricker"){
          if(fit.method=="limits"){
            # get selected sub-set of A
            #print(selectedSpec[1:2,])
            A=A[selectedSpec[,1],selectedSpec[,1]]
          }
          Acorr=sum(diag(cor(A,Aest,use="pairwise.complete.obs")))/nrow(Aest)
          oriACon=getConnectance(A)
          oriAPep=getPep(A)
          # compute number of non-zero entries correctly inferred, ignoring sign
          sens=getSensitivity(A=A,Aest=Aest)
          ppv=getPPV(A=A,Aest=Aest)
          acc=(ppv+sens)/2
          if(rand==TRUE){
            rand.correls=c()
            rand.accs=c()
            for(i in 1:randiter){
              # preserves diagonal, therefore replaced by generateRandA
              #Arand=generateA(N=nrow(A),c=round(oriACon,2))
              Arand=generateRandA(A)
              # cor.rand=sum(diag(cor(A,Arand)))/nrow(A)
              cor.rand=getCrossCorrel(A=A,Aest=Arand)
              rand.correls=c(rand.correls,cor.rand)
              sens.rand=getSensitivity(A=A,Aest=Arand)
              ppv.rand=getPPV(A=A,Aest=Arand)
              acc.rand=(sens.rand+ppv.rand)/2
              rand.accs=c(rand.accs,acc.rand)
            }
            Arandperc=mean(rand.accs)
            Arand=mean(rand.correls)
          }else{
            Arand=NA
          }
          #print(paste("mean correl A vs Aest",Acorr))
        }else{
          Acorr=NA
          oriACon=NA
          oriAPep=NA
          Arand=NA
          Arandperc=NA
          acc=NA
          sens=NA
          ppv=NA
        }
      }else{
        Acorr=NA
        oriACon=NA
        oriAPep=NA
        Arand=NA
        Arandperc=NA
        acc=NA
        sens=NA
        ppv=NA
      }
      limitsqualAcorr=c(limitsqualAcorr,Acorr)
      limitsqualmaxcorr=c(limitsqualmaxcorr,max(crosscorr1,na.rm=TRUE))
      limitsqualcrosscorAll=c(limitsqualcrosscorAll,crosscorr1All)
      limitsqualslope=c(limitsqualslope,linregslope)
      limitsqualdeltaAest=c(limitsqualdeltaAest,deltaAest)
      limitsqualautocorslope=c(limitsqualautocorslope, linregautoslope)
      limitsqualConn=c(limitsqualConn, Acon)
      limitsqualPep=c(limitsqualPep,Apep)
      limitsqualdimA=c(limitsqualdimA,AestDim)
      limitsqualsumA=c(limitsqualsumA,Asum)
      limitsqualsrandA=c(limitsqualsrandA, Arand)
      limitsqualsrandAperc=c(limitsqualsrandAperc,Arandperc)
      limitsqualAcc=c(limitsqualAcc,acc)
      limitsqualSens=c(limitsqualSens,sens)
      limitsqualPPV=c(limitsqualPPV,ppv)
      oriConn=c(oriConn,oriACon)
      oriPep=c(oriPep,oriAPep)
    }else{
      # treat gaps
      limitsqualcrosscorAll=c(limitsqualcrosscorAll,NA)
      limitsqualmaxcorr=c(limitsqualmaxcorr,NA)
      limitsqualslope=c(limitsqualslope,NA)
      limitsqualdeltaAest=c(limitsqualdeltaAest,NA)
      limitsqualPep=c(limitsqualPep,NA)
      limitsqualAcorr=c(limitsqualAcorr,NA)
      limitsqualautocorslope=c(limitsqualautocorslope,NA)
      limitsqualConn=c(limitsqualConn, NA)
      limitsqualdimA=c(limitsqualdimA,NA)
      limitsqualsumA=c(limitsqualsumA,NA)
      limitsqualsrandA=c(limitsqualsrandA, NA)
      limitsqualAcc=c(limitsqualAcc,NA)
      limitsqualsrandAperc=c(limitsqualsrandAperc,NA)
      limitsqualSens=c(limitsqualSens,NA)
      limitsqualPPV=c(limitsqualPPV,NA)
      oriConn=c(oriConn,NA)
      oriPep=c(oriPep,NA)
    }
  }

  # assemble table
  resulttable=list(limitsqualslope, limitsqualautocorslope,limitsqualcrosscorAll, limitsqualmaxcorr, limitsqualdeltaAest, limitsqualAcorr, limitsqualConn, oriConn, limitsqualPep, oriPep, limitsqualdimA, limitsqualsumA, limitsqualsrandA, limitsqualsrandAperc, limitsqualAcc,limitsqualSens,limitsqualPPV)
  names(resulttable)=c("slope","autoslope","corrAll","maxcorr","deltaAest", "Acorr","Acon","oriAcon", "Apep", "oriApep", "Adim","Asum","Arand","Arandacc","acc","sens","ppv")
  return(resulttable)
}

# compute mean cross-correlation
# between rows of A and inferred A
getCrossCorrel<-function(A,Aest){
  ccs=c()
  for(row.index in 1:nrow(A)){
    # avoid computing correlations between 2 zero rows
    if(sum(Aest[row.index,])>0 || sum(A[row.index,])>0){
      cc=cor(as.numeric(Aest[row.index,]),as.numeric(A[row.index,]), use="pairwise.complete.obs")
      ccs=c(ccs,cc)
    }
  }
  return(mean(ccs))
}

# get the percentage of correctly inferred entries in the
# original interaction matrix, ignoring sign and strength
# (percentage of non-zero entries in A that are also non-zero in Aest)
# corresponds to sensitivity
# A: original interaction matrix
# Aest: inferred interaction matrix
getSensitivity<-function(A, Aest){
  nonzero.indices=which(A!=0,arr.ind = TRUE)
  num.correct.indices=0
  for(nonzero.index in 1:nrow(nonzero.indices)){
    inferredVal=Aest[nonzero.indices[nonzero.index,1],nonzero.indices[nonzero.index,2]]
    if(!is.na(inferredVal) && inferredVal!=0){
      num.correct.indices=num.correct.indices+1
    }
  }
  # percentage of correctly inferred non-zero entries
  hundredPerc=nrow(nonzero.indices)
  onePerc=hundredPerc/100
  percent.correct.indices=num.correct.indices/onePerc
  return(percent.correct.indices)
}

# get the percentage of entries in the inferred interaction matrix
# that are correct, ignoring sign and strength
# (percentage of non-zero entries in Aest that are also non-zero in A)
# corresponds to positive predictive value
# A: original interaction matrix
# Aest: inferred interaction matrix
getPPV<-function(A, Aest){
  nonzero.indices=which(A!=0,arr.ind = TRUE)
  #print(paste("Number of links in original matrix:",nrow(nonzero.indices)))
  num.correct.indices=0
  for(nonzero.index in 1:nrow(nonzero.indices)){
    inferredVal=Aest[nonzero.indices[nonzero.index,1],nonzero.indices[nonzero.index,2]]
    if(!is.na(inferredVal) && inferredVal!=0){
      num.correct.indices=num.correct.indices+1
    }
  }
  #print(paste("Number of correctly inferred links:",num.correct.indices))
  nonzero.indices.aest=which(Aest!=0,arr.ind = TRUE)
  #print(paste("Number of inferred links:",nrow(nonzero.indices.aest)))
  # percentage of correctly inferred non-zero entries among all positive entries
  hundredPerc=nrow(nonzero.indices.aest)
  onePerc=hundredPerc/100
  percent.correct.indices=num.correct.indices/onePerc
  return(percent.correct.indices)
}
