#' Collect fitting results to community time series.
#' Different quality scores are computed for each of the three
#' fitting methods: limits, neutral and soc.
#' If the input settings are provided via the input.folder argument,
#' the parameters obtained through fitting are compared to the
#' real parameters.
#'
#' @param fit.folder the folder where fitting results are stored.
#' @param input.folder (optional) the folder where results of function generateTS are stored.
#' @param expIds the experiment identifiers of time series to be considered
#' @param fit.method the fitting method (limits, neutral or soc)
#' @param fit.type the part of the time series used (all, selected, first)
#' @return a table with goodness of fit scores

compareFit<-function(fit.folder="", input.folder="", expIds=c(), fit.method="limits", fit.type="all"){
  if(fit.folder != ""){
    if(!file.exists(fit.folder)){
      stop(paste("The fitting result folder",fit.folder,"does not exist!"))
    }
    input.timeseries.folder=paste(fit.folder,"timeseries",sep="/")
    if(!file.exists(input.timeseries.folder)){
      stop("The fitting result folder does not have a time series subfolder!")
    }
  }else{
    stop("Please provide the fitting result folder!")
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

  limitsqualslope=c()
  limitsqualmaxcorr=c()
  limitsqualdeltaAest=c()
  limitsqualAcorr=c()

  for(expId in expIds){

    print(paste("Processing identifier",expId))
    input.timeseries.name=paste(expId,"timeseries",sep="_")
    input.timeseries.expId.folder=paste(input.timeseries.folder,input.timeseries.name,sep="/")
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
      if((expId >= 25 && expId <= 30) || expId >= 44){
        subfolder=paste(expId,"ts_full_N60_skip1_tmax",sep="_")
      }else{
        if(fit.type=="all"){
          subfolder=paste(expId,"ts_N60_skip1_tmax3000",sep="_")
        }else if(fit.type=="selected"){
          # [1501:1600]
          subfolder=paste(expId,"ts_transcient_N60_skip1_tmax100",sep="_")
        }else if(fit.type=="first"){
          # [1:100]
          subfolder=paste(expId,"ts_N60_skip1_tmax100",sep="_")
        }
      }
      expId.subfolder = paste(input.timeseries.expId.folder,subfolder,sep="/")
      expId.subfolder.corrfile=paste(expId.subfolder, paste(expId, "timeseries_corr.txt", sep="_"), sep="/")
      #print(expId.subfolder.corrfile)

      # format corrfile:
      # the first line labels the columns (number of species)
      # the second the auto-correlation 1 step ahead
      # the 3rd, the auto-correlaion 2 steps ahead
      # the 4th, the auto-correlaion 3 steps ahead
      # the 5th, the auto-correlaion 4 steps ahead
      # the 6th, the auto-correlation 5 steps ahead
      # the 7th, the correlation 1 step ahead,
      # ...
      # the 11th, the correlation 5 steps ahead
      corrs=read.table(file=expId.subfolder.corrfile,header=FALSE)
      corrs=as.matrix(corrs)
      specnum=as.numeric(corrs[1,])
      #print(specnum)
      crosscorr1=as.numeric(corrs[7,])
      #print(crosscorr1)
      reg.data=data.frame(crosscorr1,specnum)
      linreg = lm(formula = crosscorr1~specnum)

      expId.subfolder.Aestfile=paste(expId.subfolder, paste(expId, "timeseries_Best.txt", sep="_"), sep="/")
      Aest=read.table(file=expId.subfolder.Aestfile,header=FALSE)
      Aest=as.matrix(Aest)
      print(paste("Number of selected species",nrow(Aest)))
      deltaAest=max(Aest)-min(Aest)

      if(input.folder != ""){
        if(Algorithm == "glv" || Algorithm == "soc" || Algorithm == "ricker"){
          # get selected sub-set of A
          expId.subfolder.selectedSpecFile=paste(expId.subfolder,paste(expId,"timeseries_name_species_kept.txt", sep="_"),sep="/")
          selectedSpec=as.matrix(read.table(file=expId.subfolder.selectedSpecFile,header=FALSE))
          #print(selectedSpec[1:2,])
          Aselect=A[selectedSpec[,1],selectedSpec[,1]]
          Acorr=sum(diag(cor(Aselect,Aest)))/nrow(Aest)
        }else{
          Acorr=0
        }
        limitsqualAcorr=c(limitsqualAcorr,Acorr)
      }

      limitsqualmaxcorr=c(limitsqualmaxcorr,max(crosscorr1,na.rm=TRUE))
      limitsqualslope=c(limitsqualslope,linreg$coefficients[2])
      limitsqualdeltaAest=c(limitsqualdeltaAest,deltaAest)
    }
  }

  # assemble table
  resulttable=list(limitsqualslope, limitsqualmaxcorr, limitsqualdeltaAest, limitsqualAcorr)
  names(resulttable)=c("slope","maxcorr","deltaAest", "Acorr")
  return(resulttable)
}
