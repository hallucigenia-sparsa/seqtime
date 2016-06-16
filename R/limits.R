#' LIMITS implementation
#
#' LIMITS is an algorithm developed by Fisher & Mehta to estimate the interaction matrix assuming a Ricker model.
#' @param x time series with taxa as rows and time points as columns
#' @param bagging.iter the number of iterations used for bagging
#' @param verbose print which taxon LIMITS is processing
#' @return the estimated interaction matrix
#' @references Fisher & Mehta (2014). Identifying Keystone Species in the Human Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression. PLoS One 9:e102451
#' @examples
#' N=20
#' A=generateA(N,c=0.1)
#' ts=ricker(N=N,A=A)
#' Aest=limits(ts,verbose=TRUE)
#' par(mfrow=c(2,1))
#' plotA(A,header="original")
#' plotA(Aest,header="estimated")
#' par(mfrow=c(1,1))
#' @export

limits<-function(x, bagging.iter=100, verbose=FALSE){
  x=t(x)
  print(paste("Time series has",ncol(x),"taxa"))
  if(verbose==TRUE){
    print("Processing first taxon.")
  }
  Aest<-t(limitscolumnwise(x, 1, r=bagging.iter))
  # loop taxa
  for(i in 2:ncol(x)){
    if(verbose==TRUE){
      print(paste("Processing taxon",i))
    }
    Aest<-cbind(Aest,t(limitscolumnwise(x, i, r=bagging.iter)))
  }
  return(t(Aest))
}

#### R is a time series with the columns corresponding to the species N, and
#### lines to different time points Ntp. If the first column of your time series
#### i.e. dim(R)=Ntp N
#### labels time, remove it before applying limits (ie, do R<-R[:,2:end] ).
#### i is the column of the interaction matrix we try to reconstruct.

#### output
#### Beval : the output of the function is the estimation of the "ith" line of the interaction matrix
#### error : mean of the errors made on evaluation of "y" using the estimated B matrix normalized by the variance of y ("one-time-step evalutaion").
#### errorF : idem as error but without the normalization by the variance.

# by Sophie de Buyl, translated from a Mathematica code provided by Charles Fisher and Pankaj Mehta.
limitscolumnwise <- function(R,i, r=100){

  listnumbkeysp<-c(); #list of number species kept. NOT MANDATORY

  # choices to be made:
  thresh <- .5 # orignially put to 5 (diminish to increase precision)

  # manipulating R to put the data in the appropriated form to apply limits:

  R[R == 0] <- 2.22507e-308  # we will have to take the log, so let's remove zero's.
  sd<-dim(R)
  N<-sd[2] #number of species
  Ntp<-sd[1]-1 #number of time points

  # comput medians column-wise
  colMedians=apply(R,2,median)
  # formulation with dependency on matrixStats
  #data<- R-t(kronecker(matrix(1,1,sd[1]),t(t(colMedians(R)))));#first N column of data matrix needed for limits
  data<- R-t(kronecker(matrix(1,1,sd[1]),colMedians))
  data<-data[1:(sd[1]-1),]
  data<-cbind(data,(log(R[2:sd[1],i])-log(R[1:(sd[1]-1),i])))

  # variable initiation
  res<-array(0,c(N,1)) # array for storing results

  listspecies<-seq(1,N) # to construct the choices of excluded/includes species

  for(k in 1:r){ # we do r times the same thing to get r estimations the interaction matrix B

    #initialize covariates to be included/excluded from regression model
    if (i!= 1 & i!=N){
      c1<-1:(i-1)
      c2<-(i+1):N
      excluded <-  c(c1,c2)
      included <- i
    }

    if(i==1){
      excluded<-2:N
      included<-i
    }
    if(i==N){
      excluded<-1:(N-1)
      included<-N
    }

    #randomly partition data into training and test sets

    trainset <- as.matrix(sample(1:Ntp,round(Ntp/2)))
    testset <- as.matrix(setdiff(1:Ntp,trainset))

    data1 <- data[trainset,]
    data2 <- data[testset,]

    test <- included
    results<- restrictedleastsquare(data1,data2,test) #perform the estimation
    errorEt <-results[[1]]
    B1t<-results[[2]]

    # loop that adds covariates(=species) as long as prediction error decreases by thresh
    xxx=1
    yyy=1

    while(xxx == 1 && yyy != N){

      yyy=yyy+1

      #loop to create list of performances for all regression models including an additional covariate

      #initial the loop
      test <- c(included,excluded[1])
      errorlist<-restrictedleastsquare(data1,data2,test)[[1]]
      #the loop
      for(kk in 2:(length(excluded)-1)){
        test = c(included,excluded[kk])
        errorlist<-c(errorlist,restrictedleastsquare(data1,data2,test)[[1]])
      }

      #sort the list so that prev[[1]] has the lowest prediction error
      ind<-which(errorlist == min(errorlist), arr.ind = TRUE)
      if(dim(as.matrix(ind))[1]>1){
        ind<-ind[[1]]
      }
      test <- c(included,excluded[ind])  #we choose the test set with the smallest error
      resulttemp <- restrictedleastsquare(data1, data2,test) #we re-run limits for that test set since we didn't store results
      errorEtemp <-resulttemp[[1]]
      B1temp<-resulttemp[[2]]
      # redefine included excluded covariates to account for including new covariate
      included <- test
      excluded <- setdiff(union(listspecies,test),intersect(listspecies,test))

      # if improvement over the current best model by thresh, include new species, otherwise exit loop
      if(is.na(errorEt)==FALSE && errorEt!=0 && (100.0*(errorEtemp - errorEt)/errorEt< -thresh)){ #we keep adding species
        errorEt <- errorEtemp# we update the error to be compared with.
        B1t <- B1temp
        testt <- test
      } else{ # we stop adding species
        xxx<-0
        listnumbkeysp<-c(listnumbkeysp,yyy)
      }
    } # end of while loop

    # store final regression coefficients in res
    B <- matrix(0,N,1)
    B[test]<-t(B1temp)

    res<-cbind(res, B)

  } #end or r loop

  # Bagging step: output median of res

  res<-res[,-1]
  # compute medians row-wise
  rowMedians=apply(res,1,median)
  Beval<-t(as.matrix(rowMedians))

  #listnumbkeysp

  return(Beval)
}

##############################################################
# restrictedleastsquare
# function need for limits.r
#
# by Sophie de Buyl, translated from a Mathematica code provided by from Charles Fisher and Pankaj Mehta.
################################################################
restrictedleastsquare <- function(inbag,outbag,test){

  # This subroutine performs the linear regression and computes the error on the test set
  #"inbag" is a matrix containing the training set. One row of the matrix looks like {x_1(t) - u_1, ..., x_N(t) -
  # u_N, ln x_i(t+1) - ln x_i(t) },
  #where u_i is the median abundance of species i. ;
  #"outbag" is a matrix containing the test set.
  #It is in the same form as inbag.;
  #"test" is a vector of integers specifying which covariates are included in the regression. I.e. if species 1,2,3 are included then test = {1,2,3}.

  #perform linear regression on training set
  temp<-dim(inbag)
  lastcol<-temp[2]
  X1 <- as.matrix(inbag[,test])
  y1 <- as.matrix(inbag[,lastcol])
  B1 <- as.matrix(MASS::ginv(X1,tol=2e-308)%*%y1)

  # calculate prediction error on test set

  X2 <- as.matrix(outbag[,test])
  y2 <- as.matrix(outbag[,lastcol])
  errorE <- mean((y2-X2%*%B1)*(y2-X2%*%B1))/var(y2)
  result <- list(errorE,B1)
  return(result)
}
