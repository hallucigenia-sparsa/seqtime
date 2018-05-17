#' @title Identify Noise Types
#' @description Identify noise types in a matrix row-wise. The noise type is the slope of the frequency versus the spectral density in log-log scale.
#' @details Row sums need to be above the given abundance threshold and the spectral density needs to scale significantly with the frequency (p-value below 0.05) in log-log scale.
#' The periodogram is computed with \code{\link{spectrum}} from the \pkg{stats} package.
#' The function returns a noisetypes object, which groups matrix row indices by noise type.
#' Note that not all rows may be assigned a noise type.
#' @param x a matrix with objects as rows and time points as columns
#' @param epsilon allowed deviation from the expected slope of 0 for white noise, -1 for pink noise and -2 for brown noise (all rows with a slope below -3 are classified as having black noise)
#' @param predef use predefined slope ranges to determine noise type from the periodogram; overrides epsilon (slope ranges are: black=[-INF,-2.25]; brown=[-1.75, -2.25), pink=(-1.75,-0.5], white=(-0.5,INF])
#' @param detrend remove a linear trend before computing the spectral density (recommended)
#' @param permut permute time points before computing noise types
#' @param pval.threshold significance threshold for periodogram powerlaw goodness of fit
#' @param smooth instead of fitting a line to the frequency-spectrum power law, fit a spline with function smooth.spline and consider the minimum of its derivative as the slope. In this case, the frequency is no longer required to scale significantly with the spectral density in log-log scale.
#' @param df smooth.spline parameter (degrees of freedom)
#' @param abund.threshold minimum sum per row
#' @param groups vector of group assignments with as many entries as x has samples, if non-empty computes frequencies and spectral densities for each group separately and compute noise type on pooled frequencies and spectral densities
#' @return S3 noisetypes object
#' @examples
#' N <- 10
#' ricker.out <- ricker(N,generateA(N),K=rep(0.01,N))
#' noise <- identifyNoisetypes(ricker.out, abund.threshold=0)
#' plotNoisetypes(noise)
#'
#' @export
identifyNoisetypes <- function(x, epsilon = 0.2, predef=FALSE, detrend=TRUE, pval.threshold = 0.05, smooth=FALSE, df=max(2,log10(ncol(x))), permut=FALSE, abund.threshold=0, groups=c()){
  if(predef==FALSE && (epsilon < 0 || epsilon > 0.5)){
    stop("Please select a value between 0 and 0.5 for epsilon.")
  }
  if(length(groups)>0){
    if(length(groups)!=ncol(x)){
      stop("Each sample should have a group assigned.")
    }
  }
  pink=c()
  brown=c()
  black=c()
  white=c()
  belowThreshold=c()
  nonclass=c()
  nonsig=c()
  if(permut==TRUE){
    indices=sample(1:ncol(x))
    x=x[,indices]
  }
  if(length(groups)>0){
    unique.groups=unique(groups)
    frequencies=list()
    spectra=list()
    firstGroup=TRUE
    # loop groups and merge group-specific frequencies and spectra
    for(group in unique.groups){
      print(paste("Processing group",group))
      group.member.indices=which(groups==group)
      print(paste("Number of group members:",length(group.member.indices)))
      x.sub=x[,group.member.indices]
      res.sub=ntTaxonLoop(x.sub,abund.threshold = abund.threshold, pval.threshold = pval.threshold, df=df, smooth=smooth, detrend=detrend, predef=predef, epsilon = epsilon, group.mode = TRUE)
      if(firstGroup){
        frequencies=res.sub$frequencies
        spectra=res.sub$spectra
        firstGroup=FALSE
      }else{
        # merge frequencies and spectra
        for(i in 1:nrow(x)){
          id=paste("taxon",i,sep="_")
          frequencies[[id]]=c(frequencies[[id]],res.sub$frequencies[[id]])
          spectra[[id]]=c(spectra[[id]],res.sub$spectra[[id]])
        }
      }
    } # end group loop
    # compute noise types for pooled frequencies and spectra
    for(i in 1:nrow(x)){
      id=paste("taxon",i,sep="_")
      # avoid empty values
      if(length(frequencies[[id]])>1){
        #print(paste("Computing noise type of taxon",id))
        nt=computeNoisetype(i=i,logfreq=frequencies[[id]],logspec=spectra[[id]], df=df, smooth=smooth, pval.threshold=pval.threshold, epsilon=epsilon, detrend=detrend, predef=predef, white=white, pink=pink, brown=brown, black=black, nonclass=nonclass, nonsig=nonsig)
        white=nt$white
        pink=nt$pink
        brown=nt$brown
        black=nt$black
        nonclass=nt$nonclass
        nonsig=nt$nonsig
      }else{
        nonclass=c(nonclass,i)
      }
    } # end compute noise types for pooled frequencies and intensities
  }else{
    res=ntTaxonLoop(x,abund.threshold = abund.threshold, df=df, pval.threshold = pval.threshold, smooth=smooth, detrend=detrend, predef=predef, epsilon = epsilon)
    white=res$white
    pink=res$pink
    brown=res$brown
    black=res$black
    nonclass=res$nonclass
    belowThreshold=res$belowThreshold
    nonsig=res$nonsig
  }
  print(paste("Number of taxa below the abundance threshold: ",length(belowThreshold)))
  print(paste("Number of taxa with non-significant power spectrum laws: ",length(nonsig)))
  print(paste("Number of taxa with non-classified power spectrum: ",length(nonclass)))
  print(paste("Number of taxa with white noise: ",length(white)))
  print(paste("Number of taxa with pink noise: ",length(pink)))
  print(paste("Number of taxa with brown noise: ",length(brown)))
  print(paste("Number of taxa with black noise: ",length(black)))
  res=noisetypes(white,pink,brown,black)
  return(res)
}

# Loop over taxa and identify the noise type of each taxon.
ntTaxonLoop<-function(x,abund.threshold=NA,  pval.threshold=NA, smooth=FALSE, df=NA, detrend=FALSE, predef=FALSE, epsilon=NA, group.mode=FALSE){
  pink=c()
  brown=c()
  black=c()
  white=c()
  belowThreshold=c()
  nonclass=c()
  nonsig=c()
  frequencies=list()
  spectra=list()
  # loop over rows of x
  for(i in 1:nrow(x)){
    id=paste("taxon",i,sep="_")
    #print(i)
    sum.taxon=sum(x[i,])
    # total abundance should be above threshold
    if(sum.taxon > abund.threshold){
      # spectrum is hidden by igraph, package name required
      out=stats::spectrum(x[i,], plot=FALSE, detrend=detrend)
      # avoid zeros and negative values
      if(length(which(out$freq<=0))==0 && length(which(out$spec<=0))==0){
        logfreq=log10(out$freq)
        logspec=log10(out$spec)
        if(group.mode){
          #print(paste("Frequencies and spectral densities have been computed for",id))
          frequencies[[id]]=logfreq
          spectra[[id]]=logspec
          #print(spectra[[id]])
        }else{
          nt=computeNoisetype(i=i,logfreq=logfreq,logspec = logspec, df=df, pval.threshold=pval.threshold, smooth=smooth, detrend=detrend, predef=predef, epsilon=epsilon, white=white,pink=pink,brown=brown,black=black,nonclass = nonclass, nonsig=nonsig)
          white=nt$white
          pink=nt$pink
          brown=nt$brown
          black=nt$black
          nonclass=nt$nonclass
          nonsig=nt$nonsig
        } # end no group mode
      }else{
        # zero or negative frequencies and powers
        nonclass=c(nonclass,i)
        if(group.mode){
          frequencies[[id]]<-c()
          spectra[[id]]<-c()
        }
      }
    } # taxon too rare
    else{
      belowThreshold=c(belowThreshold,i)
      if(group.mode){
        #print(paste("Taxon",id,"is below abundance threshold"))
        frequencies[[id]]<-c()
        spectra[[id]]<-c()
      }
    }
  } # end taxon loop
  if(group.mode){
    res=list(frequencies,spectra)
    names(res)=c("frequencies","spectra")
  }else{
    res=list(white,pink,brown,black,nonclass,nonsig,belowThreshold)
    names(res)=c("white","pink","brown","black","nonclass","nonsig","belowThreshold")
  }
  return(res)
}

# Compute the noise type given the index, the frequencies and the spectral
# densities of a taxon.
# i: current taxon index
# logfreq: logarithmic frequencies from the periodogram
# logspec: logarithmic spectral densities from the periodogram
# epsilon: allowed deviation of slope from integers
# predef: use predefined slope ranges to identify noise type
# df: degrees of freedom for smoothing
# smooth: use smoothing
# pval.threshold: threshold on the p-value
# white: vector of previously identified white noise taxa
# pink: vector of previously identified pink noise taxa
# brown: vector of previously identified brown noise taxa
# black: vector of previously identified black noise taxa
# nonclass: vector of previously encountered taxa with unclassified noise type
# nonsig: vector of previously encountered taxa where slope does not fit a line
computeNoisetype<-function(i=NA,logfreq=c(),logspec=c(), smooth=FALSE, df=NA, pval.threshold=NA, epsilon=NA, predef=FALSE, detrend=FALSE, white=c(), pink=c(), brown=c(), black=c(), nonclass=c(), nonsig=c()){
  if(smooth){
    # suggestion Beatrice
    sspline=smooth.spline(logfreq,logspec,df=df)
    deriv=predict(sspline, logfreq, deriv=1)
    slope=min(deriv$y)
    nt=noisetypeFromSlope(i=i,slope=slope,epsilon = epsilon, predef=predef, white=white, pink=pink, brown=brown, black=black, nonclass=nonclass)
    white=nt$white
    pink=nt$pink
    brown=nt$brown
    black=nt$black
    nonclass=nt$nonclass
  }else{
    reg.data=data.frame(logfreq,logspec)
    linreg = lm(formula = logspec~logfreq)
    slope=linreg$coefficients[2]
    sum=summary(linreg)
    pval=1-pf(sum$fstatistic[1], sum$fstatistic[2], sum$fstatistic[3])
    if(pval < pval.threshold){
      nt=noisetypeFromSlope(i=i,slope=slope,epsilon = epsilon, predef=predef, white=white, pink=pink, brown=brown, black=black, nonclass=nonclass)
      white=nt$white
      pink=nt$pink
      brown=nt$brown
      black=nt$black
      nonclass=nt$nonclass
    }else{
      # no line could be fitted
      nonsig=c(nonsig,i)
    }
  } # end no smooth
  res=list(white,pink,brown,black,nonclass,nonsig)
  names(res)=c("white","pink","brown","black","nonclass","nonsig")
  return(res)
}

# i: current taxon index
# slope: slope of power law in periodogram
# epsilon: allowed deviation of slope from integers
# predef: use predefined slope ranges to identify noise type
# white: vector of previously identified white noise taxa
# pink: vector of previously identified pink noise taxa
# brown: vector of previously identified brown noise taxa
# black: vector of previously identified black noise taxa
# nonclass: vector of previously encountered taxa with unclassified noise type
noisetypeFromSlope<-function(i=NA, slope=NA, epsilon=NA, predef=FALSE, white=c(), pink=c(), brown=c(), black=c(), nonclass=c()){
  if(predef==FALSE){
    # white noise:
    if(slope < epsilon && slope > -epsilon){
      white=c(white,i)
      #print("white")
      # pink noise:
    }else if(slope < -(1-epsilon) && slope > -(1+epsilon)){
      pink=c(pink,i)
      #print("pink")
      # brown noise:
    }else if(slope < -(2-epsilon) && slope > -(2+epsilon)){
      brown=c(brown,i)
      #print("brown")
      # black noise:
    }else if(slope <= -3){
      black=c(black,i)
      #print("black")
    }else{
      nonclass=c(nonclass,i)
      #print("nonclass")
    }
  }else{
    # predefined noise types (no unclassified rows possible)
    # white noise:
    if(slope > -0.5){
      white=c(white,i)
      # pink noise:
    }else if(slope <= -0.5 && slope > -1.75){
      pink=c(pink,i)
      # brown noise:
    }else if(slope <= -1.75 && slope > -2.25){
      brown=c(brown,i)
      # black noise:
    }else if(slope <= -2.25){
      black=c(black,i)
    }
  }
  res=list(white,pink,brown,black,nonclass)
  names(res)=c("white","pink","brown","black","nonclass")
  return(res)
}
