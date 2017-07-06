#' @title Get the taxonomy given OTU names and lineage information
#'
#' @description The lineage information is provided in form of a matrix,
#' which contains for each OTU identifier the taxonomic levels from column
#' 2 to 8 (domain, phylum, class, order, family, genus, species) and the entire
#' lineage in column 9. The lineage in column 9 is optional.
#' Alternatively, OTU identifiers can also be stored as column names.
#' In this case, taxonomic levels are assumed to range from column 1 to 7.
#' To differentiate between these two formats, set useRownames to false or true.
#'
#' @param selected a character vector of selected OTU identifiers
#' @param lineages a lineage table
#' @param level the taxonomic level (domain, phylum, class, order, family, genus or species)
#' @param useRownames match OTU identifiers to row names instead of the first column (in this case, taxonomic levels are assumed to range from column 1 to 7)
#' @return the taxonomy of the OTUs
#' @examples
#' data(david_stoolA_otus)
#' data(david_stool_lineages)
#' sorted=sort(apply(david_stoolA_otus,1,sum),decreasing=TRUE)[1:10]
#' getTaxonomy(names(sorted),david_stool_lineages,level="family")
#' @export

getTaxonomy<-function(selected=c(),lineages,level="class", useRownames=FALSE){
  # exact OTU identifier match
  taxa=c()
  domainIndex=2
  if(useRownames==TRUE){
    domainIndex=1
    indices=match(selected,rownames(lineages))
  }else{
    indices=match(selected,lineages[,1])
  }
  if(level=="domain"){
    taxa=lineages[indices,domainIndex]
  }else if(level=="phylum"){
    taxa=lineages[indices,(domainIndex+1)]
  }else if(level=="class"){
    taxa=lineages[indices,(domainIndex+2)]
  }
  else if(level=="order"){
    taxa=lineages[indices,(domainIndex+3)]
  }
  else if(level=="family"){
    taxa=lineages[indices,(domainIndex+4)]
  }
  else if(level=="genus"){
    taxa=lineages[indices,(domainIndex+5)]
  }
  else if(level=="species"){
    taxa=lineages[indices,(domainIndex+6)]
  }
  return(taxa)
}
