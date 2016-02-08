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
#' exactly, since the SOC, Ricker and neutral models have random components).
#' Furthermore, a plot of the time series is
#' stored in settings, for control.
#'
#' @param N species number
#' @param I individual number (only for SOC and UNTB)
#' @param tend the number of time points to generate
#' @param initAbundMode the mode in which initial abundances are generated (described in function generateAbundances)
#' @param c connectivity of interaction matrix (only for SOC, ricker and gLV)
#' @param Atype method to generate the interaction matrix (described in function generateA)
#' @param d diagonal value of the interaction matrix
#' @param PEP positive edge percentage of interaction matrix (only for SOC, ricker and gLV)
#' @param m immigration rate (only for UNTB)
#' @param algorithm how time series should be generated (can be untb, soc, ricker, glv, dm, davida or davidb)
#' @param interval sampling frequency (if 1, each sample is taken, if 2, every second sample is taken etc.)
#' @param output.folder path to result folder, if specified, input parameter and time series are exported to the output folder
#' @param output.expId the identifier of the time series generation experiment
#' @param maxIter maximal iteration number for finding a stable interaction matrix
#' @return list containing the generated interaction matrix, initial abundances (they double as immigration probabilities for soc), carrying capacities (they double as growth rates for glv), extinction probabilities, time series and settings
#' @export

generateTS<-function(N=100, I=1500, tend=100, initAbundMode=5,c=0.05, Atype="klemm",d=-1, PEP=30, m=0.02, algorithm="dm", interval=1, output.folder="", output.expId="", input.folder="", input.expId="", read.A=FALSE, read.y=FALSE, read.K=FALSE, read.extinct=FALSE, maxIter=10){

  # CONSTANTS
  sigma=0.05 # noise term in ricker
  negedge.symm=FALSE # negative interactions are not forced to be symmetric
  clique.size=5 # module size for interaction matrix generated with Klemm
  count=1000 # number of read counts per sample (for generation of Dirichlet-Multinomial time series)
  theta=0.002 # Dirichlet-Multinomial overdispersion parameter
  k=0.5 # parameter for initAbundMode=6 (evenness of geometric series, does not affect any other initAbundMode)
  tstep=1 # step size for glv
  interpolation.method="stineman" # interpolation for David data
  stability.method="ricker" # method to test stability of interaction matrix
  explosion.bound=10^8 # abundance limit above which ricker reports an explosion
  max.carrying.capacity=0.5 # maximal carrying capacity/growth rate for ricker and glv
  divisor = 3 # num positive edges/divisor are converted into negative edges during tweaking

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
    if(output.expId != ""){
      exp.settings.name=paste(output.expId,"settings",sep="_")
      exp.settings.folder=paste(settings.folder,exp.settings.name,sep="/")
      if(!file.exists(exp.settings.folder)){
        dir.create(exp.settings.folder)
      }
      exp.timeseries.name=paste(output.expId,"timeseries",sep="_")
      exp.timeseries.folder=paste(timeseries.folder,exp.timeseries.name,sep="/")
      if(!file.exists(exp.timeseries.folder)){
        dir.create(exp.timeseries.folder)
      }
    }
  } # output.folder provided

  # if non-empty, check that input folder structure exists
  if(input.folder != ""){
    if(!file.exists(input.folder)){
      stop(paste("The input folder",input.folder,"does not exist!"))
    }
    input.settings.folder=paste(input.folder,"settings",sep="/")
    if(!file.exists(input.settings.folder)){
      stop("The input folder does not have a settings subfolder!")
    }
    input.settings.name=paste(input.expId,"settings",sep="_")
    input.settings.expId.folder=paste(input.settings.folder,input.settings.name,sep="/")
    if(!file.exists(input.settings.expId.folder)){
      stop("The input settings folder does not have a subfolder for the input experiment identifier!")
    }
  }

  # generate initial abundances
  if(algorithm == "glv" || algorithm == "ricker" || algorithm == "soc" || algorithm == "dm" || algorithm == "untb"){
    if(input.folder != "" && read.y==TRUE){
      y.name=paste(input.expId,"initabund.txt",sep="_")
      input.path.y=paste(input.settings.expId.folder,y.name,sep="/")
      print(paste("Reading initial abundances from:",input.path.y,sep=" "))
      y=as.numeric(scan(input.path.y))
    }else{
      # generate initial abundances between 0 and 1
      y=generateAbundances(N=N,mode=initAbundMode, count=count, k=k)
      y=y/sum(y)
      if(output.folder != ""){
        y.name=paste(output.expId,"initabund.txt",sep="_")
        y.path=paste(exp.settings.folder,y.name,sep="/")
        write(t(y),file=y.path,ncolumns=1,sep="\t")
      }
    }
  }

  # generate carrying capacities/growth rates (needed for soc too, for tests of matrix stability)
  if(algorithm == "ricker" || algorithm == "glv" || algorithm == "soc"){
    if(input.folder != "" && read.K==TRUE){
      K.name=paste(input.expId,"capacities.txt",sep="_")
      input.path.K=paste(input.settings.expId.folder,K.name,sep="/")
      print(paste("Reading capacities from:",input.path.K,sep=" "))
      K=as.numeric(scan(input.path.K))
    }else{
      K=runif(N,min=0,max=max.carrying.capacity)
      # save carrying capacities/growth rates
      if(output.folder != ""){
        K.name=paste(output.expId,"capacities.txt",sep="_")
        K.path=paste(exp.settings.folder,K.name,sep="/")
        write(t(K),file=K.path,ncolumns=1,sep="\t")
      }
    }
  }

  # generate extinction probabilities
  if(algorithm == "soc"){
    if(input.folder != "" && read.extinct==TRUE){
      ext.name=paste(input.expId,"extinction.txt",sep="_")
      input.path.ext=paste(input.settings.expId.folder,ext.name,sep="/")
      print(paste("Reading extinction probabilities from:",input.path.ext,sep=" "))
      extinction.probab=as.numeric(scan(input.path.ext))
    }else{
      extinction.probab=runif(N,min=0,max=1)
      # save extinction probabilities
      if(output.folder != ""){
        ext.name=paste(output.expId,"extinction.txt",sep="_")
        ext.path=paste(exp.settings.folder,ext.name,sep="/")
        write(t(extinction.probab),file=ext.path,ncolumns=1,sep="\t")
      }
    }
  }

  # generate A
  iter=0
  if(algorithm == "ricker" || algorithm == "glv" || algorithm == "soc"){

    if(input.folder != "" && read.A==TRUE){
      A.name=paste(input.expId,"interactionmatrix.txt",sep="_")
      input.path.A=paste(input.settings.expId.folder,A.name,sep="/")
      print(paste("Reading extinction probabilities from:",input.path.A,sep=" "))
      A=read.table(file=input.path.A,sep="\t",header=FALSE)
      A=as.matrix(A)
    }else{
      stable = FALSE
      # try to generate a stable A
      while(stable == FALSE && iter < maxIter){
        iter = iter + 1
        A=generateA(N=N,type=Atype,c=c,d=d,pep=PEP,clique.size=clique.size,negedge.symm=negedge.symm)
        stable=testStability(A,method=stability.method, K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
        if(stable == FALSE){
          # tweak matrix
          pos.num=length(A[A>0])
          # number of tweaking trials
          trials = round(pos.num/divisor)
          for(arc.index in 1:trials){
            A=modifyA(A=A,mode="tweak")
            stable = testStability(A,method=stability.method, K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
            if(stable == TRUE){
              print(paste("It took",arc.index,"positive arcs converted into negative ones to generate a stable A."))
              break
            }
          } # end loop over positive arcs
        } # not stable
      } # iterations

      # save A
      if(stable == TRUE){
        print(paste("It took",iter,"iterations to generate a stable A."))
        # print statistics of final A
        stats=getAStats(A)

        if(output.folder != ""){
          A.name=paste(output.expId,"interactionmatrix.txt",sep="_")
          A.path=paste(exp.settings.folder,A.name,sep="/")
          write(t(A),file=A.path,ncolumns=ncol(A),sep="\t")
        }
      }else{
        stop("Could not generate a stable interaction matrix for the given parameters.")
      }
    }
  } # end generate A

  # generate time series
  print("Preparation done, running simulation...")
  if(algorithm == "ricker"){
    ts.out=ricker(N=N,A=A,sigma=sigma,y=y,tend=tend, K=K, explosion.bound=explosion.bound)
    if(ts.out == -1){
      stop("Explosion occurred!")
    }
  }else if(algorithm == "glv"){
    ts.out=glv(N=N, A=A, b=K, y=y,tstart=0, tend=tend, tstep=tstep)
  }else if(algorithm == "soc"){

    # initial abundances serve as immigration probabilities
    ts.out=soc(N=N, I=I, A=A, m.vector=y, e.vector=extinction.probab, tend=tend)
  }else if(algorithm == "untb"){
    ts.out=simuntb(N=N,y=y,m=m, tskip=0, tend=tend)
  }else if(algorithm == "dm"){
    ts.out=simCountMat(N=N,pi=y, samples=tend, counts=count, distrib="dm", theta=theta)
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
  }else{
    stop(paste("Algorithm ",algorithm," not supported",sep=""))
  }

  # sub-sample the time series
  if(interval > 1){
    ts.out=tsubsample(ts.out,interval=interval)
  }

  # assemble settings string
  settings.str=paste("Experiment identifier=",output.expId,"\n","N=",N,"\n","I=",I,"\n","steps=",tend,"\n","init abundance mode=",initAbundMode,"\n","connectance=",c,"\n","A generation method=",Atype,"\n","diagonal values of A=",d,"\n","immigration rate UNTB=",m,"\n","Positive edge percentage in A=",PEP,"\n","Algorithm=",algorithm,"\n","Sampling frequency=",interval,"\n","Iterations needed for stable A=",iter,"\n",sep="")

  if(output.folder != ""){
    # save time series
    ts.name=paste(output.expId,"timeseries.txt",sep="_")
    ts.path=paste(exp.timeseries.folder,ts.name,sep="/")
    write(t(ts.out),file=ts.path,ncolumns=ncol(ts.out),sep="\t")
    # save a plot of the time series in settings folder, for control
    ts.plot.name=paste(output.expId,"timeseriesplot.pdf",sep="_")
    ts.plot.path=paste(exp.settings.folder,ts.plot.name,sep="/")
    pdf(ts.plot.path)
    tsplot(ts.out)
    dev.off()

    # save settings
    settings.name=paste(output.expId,"settings.txt",sep="_")
    settings.path=paste(exp.settings.folder,settings.name,sep="/")
    write(settings.str,file=settings.path)
  }

  # assemble result object
  result=list(A,y,K,extinction.probab,ts.out,settings.str)
  names(result)=c("A", "y", "capacities", "extinctionprobabs","timeseries","settings")
  print("Done")

  return(result)
}


