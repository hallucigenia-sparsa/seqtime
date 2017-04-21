#' @title Self-organized instable model
#'
#' @description Generate time series with the self-organized instable (SOI) model based on Sole et al. (model B).
#'
#' @param N number of species
#' @param I number of individuals
#' @param A interaction matrix
#' @param m.vector species-specific immigration probabilities (these also determine initial abundances)
#' @param e.vector species-specific extinction probabilities
#' @param tend number of time points (i.e. the number of generations)
#' @return a matrix with species abundances as rows and time points as columns
#' @seealso \code{\link{ricker}} for the Ricker model
#' @references Sole et al. 2002 "Self-organized instability in complex ecosystems"
#' @examples
#' N=10
#' A=generateA(N,c=0.1)
#' tsplot(soi(N=N,I=1000,A=A,tend=100), header="SOI")
#' @export

soi<-function(N, I, A, m.vector=runif(N), e.vector=runif(N), tend){

  results<-generate.parameters(N, I, A, m.vector, e.vector)

  abundances_over_time<-matrix(nrow=tend,ncol=N)

  for (i in 1:tend){
    for (site in 1:length(results$sites)){
      if (results$sites[[site]] == 0){
        results<-immigration(site,results$A,results$species_list,
                             results$abundances,results$immigration_prob,
                             results$extinction_prob,results$sites,
                             results$living_species)
      } else {
        results<-extinction(site,results$A,results$species_list,
                            results$abundances,results$immigration_prob,
                            results$extinction_prob,results$sites,
                            results$living_species)
      }
      if (results$sites[[site]] != 0){
        results<-interaction(site,results$A,results$species_list,
                             results$abundances,results$immigration_prob,
                             results$extinction_prob,results$sites,
                             results$living_species)
      }
    }#closes sites

    for (j in 1:length(results$abundances)){
      abundances_over_time[i,j]<-results$abundances[[j]]
    }

  }

  return(t(abundances_over_time))

}


###############
generate.parameters<-function(N, I, A, m.vector, e.vector){
  # species_list
  species_list <- list()
  for (i in 1:N) {
    species_list[[paste("Spec_",i,sep="")]]<-0
  }

  n_spec<-c()
  for (i in 1:N) {
    n_spec<-c(n_spec, paste("Spec_", i, sep = ""))
  }
  names(m.vector)<-n_spec
  names(e.vector)<-n_spec

  # start with partially filled sited list
  # sites are filled with randomly picked species
  # according to their immigration probabilities
  sites<-c()
  n<-c()
  for (i in 1:I){
    species<-sample(1:N,1)
    if (runif(1,min=0,max=1)<m.vector[species]){
      sites[i]<-species
    } else sites[i]<-0
    n<-c(n,paste("site_",i,sep=""))
  }
  names(sites)<-n

  # generate abundance list (species vs. abundance)
  abundances<-c()
  n_spec<-c()
  for (i in 1:N) {
    abundances[i]<-sum(sites == i)
    n_spec<-c(n_spec, paste("Spec_", i, sep = ""))
  }
  names(abundances)<-n_spec

  # generate/update species_list and living_species
  living_species<-vector()
  for (i in 1:N){
    occurrance<-which(sites == i)
    species_list[[i]]<-occurrance
    if (abundances[[i]]!=0){
      living_species<-append(living_species,i)
    }
  }

  parameters <- list(A=A,species_list=species_list,
                     abundances=abundances,immigration_prob=m.vector,
                     extinction_prob=e.vector,sites=sites,
                     living_species=living_species,I=I,N=N)
}


#############
# immigration of individual from randomly chosen species to empty site
# with a probability according to immigration_prob
immigration<-function(site,A,species_list,abundances,
                      immigration_prob,extinction_prob,sites,
                      living_species){

  species<-sample(1:length(abundances),1)

  if (runif(1,min=0,max=1)<immigration_prob[[species]]){
    sites[[site]]<-species
    abundances[[species]]<-abundances[[species]]+1
    species_list[[species]]<-append(species_list[[species]], site)
    if (! is.element(species, living_species)){
      living_species<-append(living_species, species)
    }
  }

  parameters <- list(site=site,A=A,species_list=species_list,
                     abundances=abundances,immigration_prob=immigration_prob,
                     extinction_prob=extinction_prob,sites=sites,living_species=living_species)
}


#####################
# death of individual on occupied sites with a probability
# according to extinction_prob
extinction<-function(site,A,species_list,abundances,
                     immigration_prob,extinction_prob,sites,living_species){

  species<-sites[[site]]
  if (runif(1,min=0,max=1)<extinction_prob[[species]]){
    sites[[site]]<-0
    abundances[[species]]<-abundances[[species]]-1
    species_list[[species]]<-species_list[[species]][species_list[[species]] != site]
    if (abundances[[species]] == 0){
      living_species<-living_species[living_species != species]
    }
  }


  parameters <- list(site=site,A=A,species_list=species_list,
                     abundances=abundances,immigration_prob=immigration_prob,
                     extinction_prob=extinction_prob,sites=sites,
                     living_species=living_species)
}


######################
# interaction according to A
# species abundance increases or decreases
# depending on the interaction coefficients

interaction<-function(site,A,species_list,abundances,
                      immigration_prob,extinction_prob,sites,living_species){

  prob<-0

  species1<-sites[[site]]
  all_sites<-1:length(sites)

  if (length(living_species)>1){
    available_species<-living_species[living_species != species1]
    if (length(available_species) == 1){
      species2 <-available_species
    } else {species2 <-sample(available_species,1)}


    if (length(species_list[[species2]]) == 1){
      site2<-species_list[[species2]]
    } else {site2<-sample(species_list[[species2]],1)}


    ############## both interaction coefficients are positive
    ############## more positive species grows with probability according to A
    if (A[species1,species2]>0 & A[species2,species1]>0){

      site3<-sample(all_sites[all_sites != site & all_sites != site2],1)

      if (A[species1,species2] > A[species2,species1]){

        if (sites[[site3]] == 0){
          prob<-A[species1,species2]+abs(A[species2,species1])

          if (runif(1,min=0,max=1)<prob){
            sites[[site3]]<-species1
            abundances[[species1]]<-abundances[[species1]]+1
          }
        }

        else {
          species3<-sites[[site3]]

          if (A[species3,species1] < prob){

            prob<-prob-A[species3,species1]

            if (runif(1,min=0,max=1)<prob){
              sites[[site3]]<-species1
              abundances[[species1]]<-abundances[[species1]]+1
              abundances[[species3]]<-abundances[[species3]]-1
            }

          }
        }

      }


      if (A[species1,species2] < A[species2,species1]){

        if (sites[[site3]] == 0){
          prob<-A[species2,species1]+abs(A[species1,species2])

          if (runif(1,min=0,max=1)<prob){
            sites[[site3]]<-species2
            abundances[[species2]]<-abundances[[species2]]+1
          }
        }

        else {
          species3<-sites[[site3]]

          if (A[species3,species2] < prob){

            prob<-prob-A[species3,species2]

            if (runif(1,min=0,max=1)<prob){
              sites[[site3]]<-species2
              abundances[[species2]]<-abundances[[species2]]+1
              abundances[[species3]]<-abundances[[species3]]-1
            }

          }
        }

      }

    }


    ############## both interaction coefficients are negative
    ############## more negative species dies with probability according to A
    if (A[species1,species2]<0 & A[species2,species1]<0){

      if (A[species1,species2] > A[species2,species1]){

        prob<-A[species1,species2]-A[species2,species1]

        if (runif(1,min=0,max=1)<prob){
          sites[[site2]]<-species1
          abundances[[species2]]<-abundances[[species2]]-1
          abundances[[species1]]<-abundances[[species1]]+1
        }
      }

      if (A[species2,species1] > A[species1,species2]){

        prob<-A[species2,species1]-A[species1,species2]

        if (runif(1,min=0,max=1)<prob){
          sites[[site]]<-species2
          abundances[[species1]]<-abundances[[species1]]-1
          abundances[[species2]]<-abundances[[species2]]+1
        }
      }
    }



    ########## competition: one coefficient is > 0, one is < 0
    if (A[species1,species2]>0 & A[species2,species1]<0){

      prob<-A[species1,species2]-A[species2,species1]

      if (runif(1,min=0,max=1)<prob){
        sites[[site2]]<-species1
        abundances[[species1]]<-abundances[[species1]]+1
        abundances[[species2]]<-abundances[[species2]]-1
      }
    }

    if (A[species2,species1]>0 & A[species1,species2]<0){

      prob<-A[species2,species1]-A[species1,species2]

      if (runif(1,min=0,max=1)<prob){
        sites[[site]]<-species2
        abundances[[species1]]<-abundances[[species1]]-1
        abundances[[species2]]<-abundances[[species2]]+1
      }
    }


    ######### one coefficient is 0, the other is < 0
    if (A[species2,species1]==0 & A[species1,species2]<0){

      prob<--A[species1,species2]

      if (runif(1,min=0,max=1)<prob){
        sites[[site]]<-species2
        abundances[[species1]]<-abundances[[species1]]-1
        abundances[[species2]]<-abundances[[species2]]+1
      }
    }

    if (A[species1,species2]==0 & A[species2,species1]<0){

      prob<--A[species2,species1]

      if (runif(1,min=0,max=1)<prob){
        sites[[site2]]<-species1
        abundances[[species2]]<-abundances[[species2]]-1
        abundances[[species1]]<-abundances[[species1]]+1
      }
    }


    ######### one coefficient is 0, the other is > 0
    if (A[species2,species1]==0 & A[species1,species2]>0){

      site3<-sample(all_sites[all_sites != site & all_sites != site2],1)
      prob<-A[species1,species2]

      if (sites[[site3]] == 0){
        if (runif(1,min=0,max=1)<prob){
          sites[[site3]]<-species1
          abundances[[species1]]<-abundances[[species1]]+1
        }
      }
      else {
        species3<-sites[[site3]]

        if (A[species3,species1] < prob){

          prob<-prob-A[species3,species1]

          if (runif(1,min=0,max=1)<prob){
            sites[[site3]]<-species1
            abundances[[species1]]<-abundances[[species1]]+1
            abundances[[species3]]<-abundances[[species3]]-1
          }

        }
      }
    }

    if (A[species1,species2]==0 & A[species2,species1]>0){

      site3<-sample(all_sites[all_sites != site & all_sites != site2],1)
      prob<-A[species2,species1]

      if (sites[[site3]] == 0){
        if (runif(1,min=0,max=1)<prob){
          sites[[site3]]<-species2
          abundances[[species2]]<-abundances[[species2]]+1
        }
      }
      else {
        species3<-sites[[site3]]

        if (A[species3,species2] < prob){

          prob<-prob-A[species3,species2]

          if (runif(1,min=0,max=1)<prob){
            sites[[site3]]<-species2
            abundances[[species2]]<-abundances[[species2]]+1
            abundances[[species3]]<-abundances[[species3]]-1
          }

        }
      }

    }

    # generate/update species_list and living_species
    living_species<-vector()

    for (j in 1:length(species_list)){
      occurrance<-which(sites %in% j)
      species_list[[j]]<-occurrance
      if (abundances[[j]]!=0){
        living_species<-append(living_species,j)
      }
    }
  }

  parameters <- list(site=site,A=A,species_list=species_list,
                     abundances=abundances,immigration_prob=immigration_prob,
                     extinction_prob=extinction_prob,sites=sites,
                     living_species=living_species)
}

