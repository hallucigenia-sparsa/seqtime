#' @title Generate community time series
#'
#' @description Wrapper function for time series generation.
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
#' exactly, since the SOI, Ricker and neutral models have random components).
#' Furthermore, a plot of the time series and (if required by the model) the interaction
#' matrix is stored in settings, for control. Properties of the interaction matrix
#' are also stored there.
#'
#' @details On the random generation of parameter values:
#' Extinction rates are sampled from the uniform distribution
#' with values between 0 and 1. Carrying capacities/growth rates are sampled
#' from the uniform distribution with values between 0 and 0.5.
#' On the implementation of the neutral model: simHubbell is used to generate
#' neutral model time series. The 1000 first time points of the simulation are always skipped.
#' The species number of the local community is set to 1/10th of the species number in the
#' metacommunity. The local community is initialized with even abundances.
#' On the interaction matrix:
#' When the positive edge percentage is modified, symmetry is not enforced. The maximal
#' absolute interaction strength in the interaction matrix is set to 1.
#' Interaction matrix stability is tested with method "ricker", using an
#' explosion bound of 10^8. For interaction stabilizing method tweak, conversion of positive into negative
#' arcs is stopped when more than a third of the positive arcs were converted already.
#' On the stool data:
#' The stool data are processed as follows: for stool data set B, the last time point is omitted.
#' The data are rarefied to 10000 reads per sample, omitting samples with insufficient read numbers.
#' Then, using metadata on sampling days, values for missing days are interpolated with method stineman.
#' The top-abundant 100 taxa (according to their sum across samples) are retained and small negative values
#' introduced by the interpolation are set to zero.
#' Other settings:
#' For initial abundance mode 6, k is set to 0.5.
#' The simulation step size is always 1 (larger intervals are introduced after simulation).
#' Dirichet-Multinomial data are generated with a count number of 1000 per sample.
#'
#' @param N species number (for hubbell, refers to the species number in the metacommunity)
#' @param I individual number (only for SOI and hubbell)
#' @param tend the number of time points to generate
#' @param initAbundMode the mode in which initial abundances are generated (described in function generateAbundances)
#' @param c connectivity of interaction matrix (only for SOI, ricker and gLV)
#' @param clique.size parameter for interaction matrix generation (only for SOI, ricker and glv and Atype mode klemm)
#' @param sigma noise term (only for ricker)
#' @param theta Dirichlet-Multinomial overdispersion parameter (only for dm)
#' @param Atype method to generate the interaction matrix (described in function generateA)
#' @param Atweak tweaking method for instable interaction matrix, either tweak (convert positive into negative arcs until the matrix is stable or one third of the positive arcs have been converted) or schur (apply Schur decomposition, does not guarantee a stable matrix)
#' @param d diagonal value of the interaction matrix
#' @param PEP positive edge percentage of interaction matrix (only for SOI, ricker and gLV)
#' @param m immigration rate (only for hubbell)
#' @param deathrate number of individuals that are replaced at each time step (only for hubbell)
#' @param algorithm how time series should be generated (can be hubbell, soi, ricker, glv, dm, davida or davidb)
#' @param interval sampling frequency (if 1, each sample is taken, if 2, every second sample is taken etc.)
#' @param output.folder (optional) path to result folder, if specified, input parameter and time series are exported to the output folder
#' @param output.expId (optional) the identifier of the time series generation experiment
#' @param input.folder (optional) the folder from which previously generated parameters can be read
#' @param input.expId (optional) the identifier of a previous experiment whose results will be read
#' @param read.A read the interaction matrix from the given previous experiment
#' @param read.K read the carrying capacities/growth rates from the given previous experiment
#' @param read.y read the initial abundances/SOI immigration rates from the given previous experiment (for hubbell, initial abundances refer to the metacommunity)
#' @param read.extinct read the SOI extinction rates from the given previous experiment
#' @param read.ts read the previously generated time series from the given previous experiment
#' @param maxIter maximal iteration number for finding a stable interaction matrix
#' @return list containing the generated interaction matrix, initial abundances (they double as immigration probabilities for soi and metapopulation proportions for hubbell), carrying capacities (they double as growth rates for glv), extinction probabilities, time series and settings
#' @export

generateTS<-function(N=100, I=1500, tend=100, initAbundMode=5,c=0.05,clique.size=5, sigma=0.05, theta=0.002, Atype="klemm", Atweak="tweak",d=-1, PEP=30, m=0.02, deathrate=10, algorithm="dm", interval=1, output.folder="", output.expId="", input.folder="", input.expId="", read.A=FALSE, read.y=FALSE, read.K=FALSE, read.extinct=FALSE, read.ts=FALSE, maxIter=10){

  if(algorithm=="soi" || algorithm=="SOI"){
    algorithm="soc"
  }

  # CONSTANTS
  negedge.symm=FALSE # negative interactions are not forced to be symmetric
  A.max=1 # maximal interaction strength in the interaction matrix
  count=1000 # number of read counts per sample (for generation of Dirichlet-Multinomial time series)
  k=0.5 # parameter for initAbundMode=6 (evenness of geometric series, does not affect any other initAbundMode)
  tstep=1 # step size for glv
  interpolation.method="stineman" # interpolation for David data
  stability.method="ricker" # method to test stability of interaction matrix
  explosion.bound=10^8 # abundance limit above which ricker reports an explosion
  max.carrying.capacity=0.5 # maximal carrying capacity/growth rate for ricker and glv
  divisor = 3 # num positive edges/divisor are converted into negative edges during tweaking with tinker method tweak
  david.minsamplesum=10000 # minimum read number to be present in a sample from the David data
  burnin = 1000 # burnin period for hubbell, skips the transient where species number increases
  localSpecNumber=round(N/10) # local species number for hubbell is set to one tenth of the species number

  if(m < 0.2){
    burnin=5000
  }

  # empty objects
  A=matrix() # interaction matrix
  y=c() # initial abundances (they double as immigration probabilities for soc and metapopulation proportions for hubbell)
  K=c() # carrying capacities (they double as growth rates for glv)
  extinction.probab=c() # extinction probabilities
  ts.out=matrix() # time series

  # output.folder given: create output folder structure if it does not exist yet
  if(output.folder != ""){
    if(!file.exists(output.folder)){
      dir.create(output.folder)
    }
    settings.folder=file.path(output.folder,"settings")
    timeseries.folder=file.path(output.folder,"timeseries")
    if(!file.exists(settings.folder)){
      dir.create(settings.folder)
    }
    if(!file.exists(timeseries.folder)){
      dir.create(timeseries.folder)
    }

    # create experiment-specific output folders
    if(output.expId != ""){
      exp.settings.name=paste(output.expId,"settings",sep="_")
      exp.settings.folder=file.path(settings.folder,exp.settings.name)
      if(!file.exists(exp.settings.folder)){
        dir.create(exp.settings.folder)
      }
      exp.timeseries.name=paste(output.expId,"timeseries",sep="_")
      exp.timeseries.folder=file.path(timeseries.folder,exp.timeseries.name)
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
    input.settings.folder=file.path(input.folder,"settings")
    if(!file.exists(input.settings.folder)){
      stop("The input folder does not have a settings subfolder!")
    }
    input.timeseries.folder=file.path(input.folder,"timeseries")
    if(!file.exists(input.settings.folder) && read.ts==TRUE){
      stop("The input folder does not have a time series subfolder!")
    }
    input.settings.name=paste(input.expId,"settings",sep="_")
    input.settings.expId.folder=file.path(input.settings.folder,input.settings.name)
    if(!file.exists(input.settings.expId.folder)){
      stop("The input settings folder does not have a subfolder for the input experiment identifier!")
    }
    input.timeseries.name=paste(input.expId,"timeseries",sep="_")
    input.timeseries.expId.folder=file.path(input.timeseries.folder,input.timeseries.name)
    if(!file.exists(input.timeseries.expId.folder) && read.ts==TRUE){
      stop("The input time series folder does not have a subfolder for the input experiment identifier!")
    }
  }

  # generate initial abundances
  if(algorithm == "glv" || algorithm == "ricker" || algorithm == "soc" || algorithm == "dm" || algorithm == "hubbell"){
    if(read.ts == FALSE){
      if(input.folder != "" && read.y==TRUE){
        y.name=paste(input.expId,"initabund.txt",sep="_")
        input.path.y=file.path(input.settings.expId.folder,y.name)
        print(paste("Reading initial abundances from:",input.path.y,sep=" "))
        y=as.numeric(scan(input.path.y))
      }else{
        # generate initial abundances between 0 and 1
        y=generateAbundances(N=N,mode=initAbundMode, count=count, k=k, probabs=TRUE)
        if(output.folder != ""){
          y.name=paste(output.expId,"initabund.txt",sep="_")
          y.path=file.path(exp.settings.folder,y.name)
          write(t(y),file=y.path,ncolumns=1,sep="\t")
        }
      }
    } # time series are not read in
  }

  # generate carrying capacities/growth rates (needed for soc too, for tests of matrix stability)
  if(algorithm == "ricker" || algorithm == "glv" || algorithm == "soc"){
    if(read.ts == FALSE){
      if(input.folder != "" && read.K==TRUE){
        K.name=paste(input.expId,"capacities.txt",sep="_")
        input.path.K=file.path(input.settings.expId.folder,K.name)
        print(paste("Reading capacities from:",input.path.K,sep=" "))
        K=as.numeric(scan(input.path.K))
      }else{
        K=runif(N,min=0,max=max.carrying.capacity)
        # save carrying capacities/growth rates
        if(output.folder != ""){
          K.name=paste(output.expId,"capacities.txt",sep="_")
          K.path=file.path(exp.settings.folder,K.name)
          write(t(K),file=K.path,ncolumns=1,sep="\t")
        }
      }
    } # time series are not read in
  }

  # generate extinction probabilities
  if(algorithm == "soc" && read.ts == FALSE){
    if(input.folder != "" && read.extinct==TRUE){
      ext.name=paste(input.expId,"extinction.txt",sep="_")
      input.path.ext=file.path(input.settings.expId.folder,ext.name)
      print(paste("Reading extinction probabilities from:",input.path.ext,sep=" "))
      extinction.probab=as.numeric(scan(input.path.ext))
    }else{
      extinction.probab=runif(N,min=0,max=1)
      # save extinction probabilities
      if(output.folder != ""){
        ext.name=paste(output.expId,"extinction.txt",sep="_")
        ext.path=file.path(exp.settings.folder,ext.name)
        write(t(extinction.probab),file=ext.path,ncolumns=1,sep="\t")
      }
    }
  }

  # generate A
  iter=0
  if(algorithm == "ricker" || algorithm == "glv" || algorithm == "soc"){

    if(input.folder != "" && read.A==TRUE){
      if(read.ts == FALSE){
        A.name=paste(input.expId,"interactionmatrix.txt",sep="_")
        input.path.A=file.path(input.settings.expId.folder,A.name)
        print(paste("Reading interaction matrix from:",input.path.A,sep=" "))
        A=read.table(file=input.path.A,sep="\t",header=FALSE)
        A=as.matrix(A)
      }
    }else{
      stable = FALSE
      # try to generate a stable A
      while(stable == FALSE && iter < maxIter){
        print(paste("Generating A, iter",iter))
        iter = iter + 1
        A=generateA(N=N,type=Atype,c=c,d=d,pep=PEP,clique.size=clique.size,negedge.symm=negedge.symm, max.strength=A.max)
        stable=testStability(A,method=stability.method, K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
        if(stable == FALSE){
          if(Atweak == "tweak"){
            # tweak matrix until it is stable or too many arcs have become positive
            pos.num=length(A[A>0])
            # number of tweaking trials
            trials = round(pos.num/divisor)
            for(arc.index in 1:trials){
              A=modifyA(A=A,mode=Atweak)
              stable = testStability(A,method=stability.method, K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
              if(stable == TRUE){
                print(paste("It took",arc.index,"positive arcs converted into negative ones to generate a stable A."))
                break
              }
            } # end loop over positive arcs
          }else if(Atweak == "schur"){
            # tweak matrix using Schur decomposition
            print("Applying Schur decomposition...")
            A=modifyA(A=A,mode=Atweak)
            stable = testStability(A,method=stability.method, K=K, y=y, sigma=sigma, explosion.bound=explosion.bound)
            if(stable == TRUE){
              print(paste("Schur decomposition made interaction matrix stable."))
              break
            }
          }else{
            stop(paste("Interaction matrix modification method",Atweak,"is not supported"))
          }
        } # not stable
      } # iterations

      # save A, its image and its properties
      if(stable == TRUE){
        print(paste("It took",iter,"iterations to generate a stable A."))
        # print statistics of final A
        stats=getAStats(A)

        if(output.folder != ""){
          A.name=paste(output.expId,"interactionmatrix.txt",sep="_")
          A.path=file.path(exp.settings.folder,A.name)
          write(t(A),file=A.path,ncolumns=ncol(A),sep="\t")
          Aplot.name=paste(output.expId,"interactionmatrix.pdf",sep="_")
          Aplot.path=file.path(exp.settings.folder,Aplot.name)
          pdf(Aplot.path)
          plotA(A,header=paste("interaction matrix experiment",output.expId))
          dev.off()
          Aprop.name=paste(output.expId,"interactionmatrix_properties.txt",sep="_")
          Aprop.path=file.path(exp.settings.folder,Aprop.name)
          Astats=getAStats(A)
          Aconnectance=getConnectance(A)
          APep=round(getPep(A),2)
          ANumArcs=length(A[A!=0])
          maxA=max(A)
          minA=min(A)
          Aprops.str=paste("Maximum value=",maxA,"\n","Minimum value=",minA,"\n","Arc number=",ANumArcs,"\n","PEP=",APep,"\n","Connectance=",Aconnectance,"\n","Number of interactions=",Astats$nbinteractions,"\n","Mean interaction strength=",round(Astats$meanstrength,4),"\n","Number of mutualisms=",Astats$nbmut,"\n","Number of competitions=",Astats$nbcomp,"\n","Number of exploitations=",Astats$nbexp,"\n","Number of commensalisms=",Astats$nbcom,"\n","Number of amensalisms=",Astats$nbam,sep="")
          write(Aprops.str,file=Aprop.path)
        }
      }else{
        stop("Could not generate a stable interaction matrix for the given parameters.")
      }
    }
  } # end generate A

  # generate time series
  if(input.folder != "" && read.ts==TRUE){
    ts.name=paste(input.expId,"timeseries.txt",sep="_")
    input.path.ts=file.path(input.timeseries.expId.folder,ts.name)
    print(paste("Reading time series from:",input.path.ts,sep=" "))
    ts.out=read.table(file=input.path.ts,sep="\t",header=FALSE)
    ts.out=as.matrix(ts.out)
  }else{
    print("Preparation done, running simulation...")
    if(algorithm == "ricker"){
      ts.out=ricker(N=N,A=A,sigma=sigma,y=y,tend=tend, K=K, explosion.bound=explosion.bound)
      if(ts.out[[1]] == -1){
        stop("Explosion occurred!")
      }
    }else if(algorithm == "glv"){
      ts.out=glv(N=N, A=A, b=K, y=y,tstart=0, tend=tend, tstep=tstep)
    }else if(algorithm == "soc"){
      # initial abundances serve as immigration probabilities
      ts.out=soi(N=N, I=I, A=A, m.vector=y, e.vector=extinction.probab, tend=tend)
    }else if(algorithm == "hubbell"){
      # after reading, for numeric reasons, may not always sum exactly to one
      y=y/sum(y)
      # local abundance distribution is even (y=rep(1/N,N))
      ts.out=simHubbell(N=localSpecNumber,M=N, I=I, d=deathrate, m.vector=y,m=m, tskip=burnin, tend=(tend+burnin))
    }else if(algorithm == "dm"){
      ts.out=simCountMat(N=N,pi=y, samples=tend, counts=count, distrib="dm", theta=theta)
    }else if(algorithm == "davida" || algorithm == "davidb"){
      if(algorithm == "davida"){
        david_stoolA_otus=matrix()
        david_stoolA_metadata=matrix()
        load(system.file("data/david_stoolA_otus.rda", package = "seqtime"))
        load(system.file("data/david_stoolA_metadata.rda", package = "seqtime"))
        stool=david_stoolA_otus
        metadata=david_stoolA_metadata
      }else if(algorithm == "davidb"){
        david_stoolB_otus=matrix()
        david_stoolB_metadata=matrix()
        load(system.file(file.path("data","david_stoolB_otus.rda"), package = "seqtime"))
        load(system.file(file.path("data","david_stoolB_metadata.rda"), package = "seqtime"))
        stool=david_stoolB_otus
        metadata=david_stoolB_metadata
        # omit the last sample in David data stool B, because there is a huge sampling gap of 66 days
        stool=stool[,1:(ncol(stool)-1)]
        metadata=metadata[,1:(ncol(metadata)-1)]
      }
      # rarefy stool data
      rarefyRes=rarefyFilter(stool,min=david.minsamplesum)
      stool=rarefyRes$rar
      # discard days with read count below minsamplesum
      days=metadata[1,rarefyRes$colindices]
      # make data equidistant
      stool.interp=interpolate(stool,time.vector=days,interval=1,method=interpolation.method)
      # select top N abundant taxa
      if(nrow(stool.interp) >= N){
        sorted=sort(apply(stool.interp,1,sum),decreasing=TRUE,index.return=TRUE)
        stool.interp=stool.interp[sorted$ix[1:N],]
      }
      # treat negative values introduced through interpolation
      negValues=stool.interp[stool.interp<0]
      meanNegVal=mean(negValues)
      rangeNegVal=range(negValues)
      print(paste("The interpolation introduced ",length(negValues)," negative values with a mean value of ",meanNegVal," a minimum of ",rangeNegVal[1]," and a maximum of ",rangeNegVal[2],". These are now set to zero.",sep=""))
      stool.interp[stool.interp<0]=0
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
  }

  # sub-sample the time series
  if(interval > 1){
    ts.out=tsubsample(ts.out,interval=interval)
  }

  # assemble settings string
  if(input.expId == ""){
    input.expId == "NA"
  }
  settings.str=paste("Output_experiment_identifier=",output.expId,"\n","Input_experiment_identifier=",input.expId,"\n","Interaction_matrix_read_from_input=",read.A,"\n","Initial_abundances_or_immigration_rates_read_from_input=",read.y,"\n","Growth_rates_or_capacities_read_from_input=",read.K,"\n","Extinction_rates_read_from_input=",read.extinct,"\n","Timeseries_read_from_input=",read.ts,"\n","N=",N,"\n","I=",I,"\n","steps=",tend,"\n","init_abundance_mode=",initAbundMode,"\n","sigma=",sigma,"\n","theta=",theta,"\n","clique_size_in_interaction_matrix=",clique.size,"\n","connectance=",c,"\n","Interaction_matrix_generation_method=\"",Atype,"\"\n","Interaction_matrix_tweaking_method=\"",Atweak,"\"\n","diagonal_values_of_interaction_matrix=",d,"\n","immigration_rate_Hubbell=",m,"\n","deathrate_Hubbell=",deathrate,"\n","Requested_positive_edge_percentage_of_interaction_matrix=",PEP,"\n","Algorithm=\"",algorithm,"\"\n","Sampling_frequency=",interval,"\n","Iterations_needed_for_stable_interaction_matrix=",iter,"\n",sep="")

  if(output.folder != ""){
    # save time series
    ts.name=paste(output.expId,"timeseries.txt",sep="_")
    ts.path=file.path(exp.timeseries.folder,ts.name)
    write(t(ts.out),file=ts.path,ncolumns=ncol(ts.out),sep="\t")
    # save a plot of the time series in settings folder, for control
    ts.plot.name=paste(output.expId,"timeseriesplot.pdf",sep="_")
    ts.plot.path=file.path(exp.settings.folder,ts.plot.name)
    pdf(ts.plot.path)
    tsplot(ts.out,header=paste("experiment",output.expId))
    dev.off()

    # save settings
    settings.name=paste(output.expId,"settings.txt",sep="_")
    settings.path=file.path(exp.settings.folder,settings.name)
    write(settings.str,file=settings.path)
  }

  # assemble result object
  result=list(A,y,K,extinction.probab,ts.out,settings.str)
  names(result)=c("A", "y", "capacities", "extinctionprobabs","timeseries","settings")
  print("Done")

  return(result)
}


