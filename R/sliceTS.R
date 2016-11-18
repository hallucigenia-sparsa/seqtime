#' Store slices of previously generated time series in the input folder
#' and store them in the output folder.
#'
#' @param slices a matrix with rows and 2 columns, where each row codes the time series identifier, the first column the start and the second the end time point, NA means until end of the time series
#' @param input.folder location of time series and settings sub folders
#' @param output.folder folder in which all slices go (no sub-folders)
#' @param expIds set of experiment identifiers to process

sliceTS<-function(slices, input.folder="", output.folder="", expIds=c()){
  if(input.folder != ""){
    if(!file.exists(input.folder)){
      stop(paste("The input folder",input.folder,"does not exist!"))
    }
    input.timeseries.folder=paste(input.folder,"timeseries",sep="/")
    if(!file.exists(input.timeseries.folder)){
      stop("The input folder does not have a time series subfolder!")
    }
  }else{
    stop("Please provide the input folder!")
  }

  if(output.folder != ""){
    if(!file.exists(output.folder)){
      dir.create(output.folder)
    }
  }else{
    stop("Please provide the output folder!")
  }

  for(expId in expIds){
    print(paste("Processing identifier",expId))

    input.timeseries.name=paste(expId,"timeseries",sep="_")
    input.timeseries.expId.folder=paste(input.timeseries.folder,input.timeseries.name,sep="/")
    if(!file.exists(input.timeseries.expId.folder)){
      stop("The input time series folder does not have a subfolder for the input experiment identifier!")
    }

    # read time series file
    ts.name=paste(expId,"timeseries.txt",sep="_")
    input.path.ts=paste(input.timeseries.expId.folder,ts.name,sep="/")
    print(paste("Reading time series from:",input.path.ts,sep=" "))
    ts=read.table(file=input.path.ts,sep="\t",header=FALSE)
    ts=as.matrix(ts)
    N=nrow(ts)

    startSlice=slices[expId,1]
    endSlice=slices[expId,2]
    if(is.na(endSlice)){
      endSlice=ncol(ts)
    }

    slicedTS=ts[,startSlice:endSlice]
    print(paste("Length of the slice",ncol(slicedTS)))

    # save time series
    ts.name=paste(expId,"sliced_timeseries.txt",sep="_")
    ts.path=paste(output.folder,ts.name,sep="/")
    write(t(slicedTS),file=ts.path,ncolumns=ncol(slicedTS),sep="\t")

  } # end loop expIds

}
