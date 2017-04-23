#' @title Get the taxonomy given OTU names and lineage information
#'
#' @description The lineage information is provided in form of a matrix,
#' which contains for each OTU identifier the taxonomic levels from column
#' 2 to 8 (domain, phylum, class, order, family, genus, species) and the entire
#' lineage in column 9.
#'
#' @param selected a character vector of selected OTU identifiers
#' @param lineages a lineage table
#' @param level the taxonomic level (domain, phylum, class, order, family, genus or species)
#' @return the taxonomy of the OTUs
#' @examples
#' data(david_stoolA_otus)
#' data(david_stool_lineages)
#' sorted=sort(apply(david_stoolA_otus,1,sum),decreasing=TRUE)[1:10]
#' getTaxonomy(names(sorted),david_stool_lineages,level="family")
#' @export

getTaxonomy<-function(selected=c(),lineages,level="class"){
  # exact OTU identifier match
  indices=match(selected,lineages[,1])
  taxa=c()
  if(level=="domain"){
    taxa=lineages[indices,2]
  }else if(level=="phylum"){
    taxa=lineages[indices,3]
  }else if(level=="class"){
    taxa=lineages[indices,4]
  }
  else if(level=="order"){
    taxa=lineages[indices,5]
  }
  else if(level=="family"){
    taxa=lineages[indices,6]
  }
  else if(level=="genus"){
    taxa=lineages[indices,7]
  }
  else if(level=="species"){
    taxa=lineages[indices,8]
  }
  return(taxa)
}
