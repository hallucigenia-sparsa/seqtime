#' Generate community time series
#'
#' Wrapper function for time series generation.
#' The result folder has a specific structure that is
#' created if it does not yet exist. It has two
#' sub-folders, one named timeseries and the other
#' settings. In both sub-folders, there are folders
#' having as name the experiment identifier followed by
#' timeseries or settings. The interaction matrix,
#' immigration and extinction rates, carrying capacities
#' and other randomly generated input are also
#' stored in settings. The idea is that the time
#' series can be re-created from the information
#' stored in the settings folder (though not
#' exactly, since all 3 models have random components).
#'
#' @param N species number
#' @param I individual number (only for SOC and UNTB)
#' @param tend the number of time points to generate
#' @param initAbundMode the mode in which initial abundances are generated (described in generateAbundances.R)
#' @param c connectivity of interaction matrix (only for SOC, ricker and gLV)
#' @param PEP positive edge percentage of interaction matrix (only for SOC, ricker and gLV)
#' @param m immigration rate (only for UNTB)
#' @param algorithm how time series should be generated (can be untb, soc, ricker, glv, dm, davida or davidb)
#' @param interval sampling frequency (if 1, each sample is taken, if 2, every second sample is taken etc.)
#' @param output.folder path to result folder, if specified, input parameter and time series are exported to the output folder
#' @param expId the identifier of the time series generation experiment
#' @param maxIter maximal iteration number for finding a stable interaction matrix (if negative or 0, no constraint is set on the iteration number)
#' @return list containing the generated interaction matrix, initial abundances (they double as immigration probabilities for soc), carrying capacities (they double as growth rates for glv), extinction probabilities, time series and settings
#' @export

generateTS<-function(N=100, I=1500, tend=2000, initAbundMode=5,c=0.02, Atype="klemm",d=-1, PEP=30, m=0.05, algorithm="", interval=1, output.folder="", expId="", maxIter=10){

  # CONSTANTS
  sigma=0.05 # noise term in ricker
  negedge.symm=FALSE # negative interactions are not forced to be symmetric
  clique.size=5 # interaction matrix generated with Klemm
  count=10000 # number of read counts per sample (for generation of initial abundances and Dirichlet-Multinomial)
  theta=0.002 # Dirichlet-Multinomial overdispersion parameter
  k=0.5 # parameter for initial abundance generation mode 6
  tstep=1 # step size for glv
  interpolation.method="stineman"

  # empty objects
  A=matrix() # interaction matrix
  y=c() # initial abundances (they double as immigration probabilities for soc)
  K=c() # carrying capacities (they double as growth rates for glv)
  extinction.probab=c() # extinction probabilities
  ts.out=matrix() # time series

  # output.folder given: create output folder structure if it does not exist yet
  if(output.folder != ""){
    if(!file.exists(output.folder)){
      dir.create(output.folder)
    }
    settings.folder=paste(output.folder,"settings",sep="/")
    timeseries.folder=paste(output.folder,"timeseries",sep="/")
    if(!file.exists(settings.folder)){
      dir.create(settings.folder)
    }
    if(!file.exists(timeseries.folder)){
      dir.create(timeseries.folder)
    }

    # create experiment-specific output folders
    if(expId != ""){
      exp.settings.name=paste(expId,"settings",sep="_")
      exp.settings.folder=paste(settings.folder,exp.settings.name,sep="/")
      if(!file.exists(exp.settings.folder)){
        dir.create(exp.settings.folder)
      }
      exp.timeseries.name=paste(expId,"timeseries",sep="_")
      exp.timeseries.folder=paste(timeseries.folder,exp.timeseries.name,sep="/")
      if(!file.exists(exp.timeseries.folder)){
        dir.create(exp.timeseries.folder)
      }
    }
  } # output.folder provided

  # generate A
  iter=1
  if(algorithm == "ricker" || algorithm == "glv" || algorithm == "soc"){
      # clicque.size defaults to 5, negative edges are by default not symmetric
      A=generateA(N=N,type=Atype,c=c,d=d,pep=PEP)
      # re-generate A until it is stable
      while(testStability(A,method="eigen")==FALSE){
        A=generateA(N=N,type=Atype,c=c,d=d,pep=PEP,clique.size=clique.size,negedge.symm=negedge.symm)
        if(maxIter > 0 && iter >= maxIter){
          stop("Cannot generate a stable interaction matrix for the given parameters")
        }
        iter = iter+1
      }
      print(paste("It took",iter,"iterations to generate a stable A."))
      # save A
      if(output.folder != ""){
        A.name=paste(expId,"interactionmatrix.txt",sep="_")
        A.path=paste(exp.settings.folder,A.name,sep="/")
        write(t(A),file=A.path,ncolumns=ncol(A),sep="\t")
      }
  }

  # generate initial abundances
  if(algorithm == "glv" || algorithm == "ricker" || algorithm == "soc" || algorithm == "dm"){
    # generate initial abundances
    y=generateAbundances(N=N,mode=initAbundMode, count=count, k=k)
    if(output.folder != ""){
      y.name=paste(expId,"initabund.txt",sep="_")
      y.path=paste(exp.settings.folder,y.name,sep="/")
      write(t(y),file=y.path,ncolumns=1,sep="\t")
    }
  }

  # generate carrying capacities/growth rates
  if(algorithm == "ricker" || algorithm == "glv"){
    K=runif(N,min=0,max=1)
    # save carrying capacities/growth rates
    if(output.folder != ""){
      K.name=paste(expId,"capacities.txt",sep="_")
      K.path=paste(exp.settings.folder,K.name,sep="/")
      write(t(K),file=K.path,ncolumns=1,sep="\t")
    }
  }

  # generate time series
  if(algorithm == "ricker"){
    ts.out=ricker(N=N,A=A,sigma=sigma,y=y,tend=tend, K=K)
  }else if(algorithm == "glv"){
    ts.out=glv(N=N, A=A, b=K, y=y,tstart=0, tend=tend, tstep=tstep)
  }else if(algorithm == "soc"){
    extinction.probab=runif(N,min=0,max=1)
    # initial abundances serve as immigration probabilities
    immigration.probab=y
    if(sum(immigration.probab) > 1){
      immigration.probab=immigration.probab/sum(immigration.probab)
    }
    ts.out=soc(N=N, I=I, A=A, m.vector=immigration.probab, e.vector=extinction.probab, tend=tend)
    # save extinction probabilities
    if(output.folder != ""){
      ext.name=paste(expId,"extinction.txt",sep="_")
      ext.path=paste(exp.settings.folder,ext.name,sep="/")
      write(t(extinction.probab),file=ext.path,ncolumns=1,sep="\t")
    }
  }else if(algorithm == "untb"){
    ts.out=simuntb(N=N,y=y,m=m, tskip=0, tend=tend)
  }else if(algorithm == "dm"){
    pi=y
    if(sum(pi) > 1){
      pi=pi/sum(pi)
    }
    ts.out=simCountMat(N=N,pi=pi, samples=tend, counts=count, distrib="dm", theta=theta)
  }else if(algorithm == "davida" || algorithm == "davidb"){
    if(algorithm == "davida"){
      data("david_stoolA_otus")
      data("david_stoolA_metadata")
      stool=david_stoolA_otus
      metadata=david_stoolA_metadata
    }else if(algorithm == "davidb"){
      data("david_stoolB_otus")
      data("david_stoolB_metadata")
      stool=david_stoolB_otus
      metadata=david_stoolB_metadata
    }
    # make data equidistant
    days=metadata[1,]
    stool.interp=interpolate(stool,time.vector=days,interval=interval,method=interpolation.method)
    # starting with the first time point, select as many samples as time points requested
    if(tend <= ncol(stool.interp)){
      ts.out=stool.interp[,1:tend]
    }else{
      # if more time points are requested than present, return the entire time series
      ts.out=stool.interp
    }
  }

  # sub-sample the time series
  if(interval > 1){
    ts.out=tsubsample(ts.out,interval=interval)
  }

  # assemble settings string
  settings.str=paste("Experiment identifier=",expId,"\n","N=",N,"\n","I=",I,"\n","steps=",tend,"\n","init abundance mode=",initAbundMode,"\n","connectance=",c,"\n","A generation method=",Atype,"\n","diagonal values of A=",d,"\n","immigration rate UNTB=",m,"\n","Positive edge percentage in A=",PEP,"\n","Algorithm=",algorithm,"\n","Sampling frequency=",interval,"\n","Iterations needed for stable A=",iter,"\n",sep="")

  if(output.folder != ""){
    # save time series
    ts.name=paste(expId,"timeseries.txt",sep="_")
    ts.path=paste(exp.timeseries.folder,ts.name,sep="/")
    write(t(ts.out),file=ts.path,ncolumns=ncol(ts.out),sep="\t")

    # save settings
    settings.name=paste(expId,"settings.txt",sep="_")
    settings.path=paste(exp.settings.folder,settings.name,sep="/")
    write(settings.str,file=settings.path)
  }

  # assemble result object
  result=list(A,y,K,extinction.probab,ts.out,settings.str)
  names(result)=c("A", "y", "capacities", "extinctionprobabs","timeseries","settings")
  print("Done")

  return(result)
}



