#' @title Compute edge statistics for an interaction matrix with taxonomic information
#'
#' @description The function expects an interaction matrix with taxon names given as row and column names
#' and their assignment to a higher taxonomic level given in the first column (can be generated with assignTaxonLevelsToA).
#' Either the edge number or the p-value of the edge number given the number of taxa in both higher-level taxon groups is reported.
#' The p-value is computed by randomizing the interaction matrix while preserving the total number of its entries.
#' Entries (i,j) and (j,i) in the interaction matrix are counted as one interaction by default.
#'
#' @param taxonA an interaction matrix with taxonomy (first column contains higher-level taxa)
#' @param count.entries instead of interactions, count links (so entries (i,j) and (j,i) are not combined, but counted separately)
#' @param count.loops count link of a node to itself
#' @param report.pval if true, entries in the resulting interaction matrix are p-values of edge number between taxonomic groups, else edge numbers
#' @param pos.only only consider positive entries in taxonA
#' @param neg.only only consider negative entries in taxonA
#' @param iters number of randomizations for report.pval
#' @param p.threshold only report higher-level links with p-values below given threshold
#' @param sig convert p-values into significances
#'
#' @return an interaction matrix with as many rows and columns as there are distinct higher-level taxa in the first column of the input
#'
#' @examples
#' data("david_stool_lineages")
#' # make a random OTU interaction matrix as an example
#' N=30
#' A=generateA(N=N,c=0.1)
#' randomOTUIndices=sample(1:nrow(david_stool_lineages))[1:N]
#' rownames(A)=david_stool_lineages[randomOTUIndices,1]
#' colnames(A)=rownames(A)
#' A.otu=assignTaxonLevelsToA(A,lineages=david_stool_lineages, taxon.level="genus",no.merge=TRUE)
#' A.phylum=analyseTaxaInA(A.otu,count.loops=TRUE, count.entries=TRUE)
#' @export

analyseTaxaInA<-function(taxonA, count.entries=FALSE, count.loops=FALSE, report.pval=FALSE, pos.only=FALSE, neg.only=FALSE, iters=100, p.threshold=0.05, sig=FALSE){
  taxon.assignment=taxonA[,1]
  A=matrix(as.numeric(taxonA[,2:ncol(taxonA)]),nrow=nrow(taxonA),ncol=nrow(taxonA))
  if(count.loops==FALSE){
    diag(A)=0 # set self-arcs to zero
  }
  higherTaxa=unique(taxon.assignment)
  A.up=matrix(NA,nrow=length(higherTaxa),ncol=length(higherTaxa))

  # only count negative or positive interactions
  if(pos.only==TRUE){
    A[A<0]=0
  }else if(neg.only==TRUE){
    A[A>0]=0
  }

  for(i in 1:length(higherTaxa)){
    for(j in 1:i){
      membersGroup1=which(taxon.assignment==higherTaxa[i])
      membersGroup2=which(taxon.assignment==higherTaxa[j])
      A.up[i,j]=countInteractions(A,group1.indices = membersGroup1, group2.indices = membersGroup2, count.entries = count.entries, count.loops = count.loops)
      A.up[j,i]=A.up[i,j]
    }  # inner higher taxon loop
  }  # outer higher taxon loop

  if(report.pval==TRUE){

    # stores number of times that randomized matrix generated more interactions than original matrix
    tempMat=matrix(1,nrow=nrow(A.up),ncol=ncol(A.up))

    for(iter in 1:iters){
      #print(paste("iter",iter))
      randA=generateRandA(A, count.loops=count.loops)
      for(i in 1:length(higherTaxa)){
        # higher-level taxon matrix is symmetric
        for(j in 1:i){
          membersGroup1=which(taxon.assignment==higherTaxa[i])
          membersGroup2=which(taxon.assignment==higherTaxa[j])
          rand.num=countInteractions(randA,group1.indices = membersGroup1, group2.indices = membersGroup2, count.entries = count.entries, count.loops = count.loops)
          if(rand.num >= A.up[i,j]){
            tempMat[i,j]=tempMat[i,j]+1
            tempMat[j,i]=tempMat[i,j]
          }
        } # inner higher taxon loop
      } # outer higher taxon loop
    } # loop iterations

    # compute p-values parameter-free
    A.up=tempMat/(iters+1)  # (rand.num.greater+1)/(iters+1)

    A.up[A.up>=p.threshold]=1

    if(sig==TRUE){
      A.up=-1*log10(A.up)
    }

  } # end compute p-values

  rownames(A.up)=higherTaxa
  colnames(A.up)=rownames(A.up)
  return(A.up)
}

# generate a random interaction
# matrix that preserves the connectance
# of the original
generateRandA<-function(A, count.loops=FALSE){
  N=nrow(A)
  # do not count diagonal
  if(count.loops==FALSE){
    diag(A)=0
  }
  targetEdgeNum=length(A[A!=0])
  Arand=matrix(0,nrow=nrow(A),ncol=ncol(A))
  for(edge.index in 1:targetEdgeNum){
    rand.x=sample(1:N)[1]
    rand.y=sample(1:N)[2]
    Arand[rand.x,rand.y]=runif(1,-0.5,0.5)
  }
  return(Arand)
}


# given group memberships, count between-group interactions
countInteractions<-function(A, group1.indices=c(), group2.indices=c(), count.entries=FALSE, count.loops=FALSE){
  num.interactions=0
  end=NA

  for(i in 1:ncol(A)){
    if(count.entries==TRUE){
      end=ncol(A)
    }else{
      end=i
    }
    # skip upper triangle
    for(j in 1:end){
      # only count loops if requested
      if((i != j) || count.loops==TRUE){
        # group indices
        if((i %in% group1.indices && j %in% group2.indices) || ((j %in% group1.indices && i %in% group2.indices))){
          # collect non-zero interaction strengths
          if(count.entries==TRUE){
            if(A[i,j]!=0){
              num.interactions=num.interactions+1
            }
          }else{
            if((A[i,j] != 0) || (A[j,i] != 0)){
              num.interactions=num.interactions+1
            }
          }
        }
      } # avoid diagonal
    } # j loop
  } # i loop
  return(num.interactions)
}
