#' David et al. stool OTU lineages
#'
#' The lineages for OTUs in David data stool subject A and B.
#' The first column holds the OTU, column 2 to 8 hold the domain, phylum,
#' class, order, family, genus and species level classification and column 9 holds the lineage
#' ending with the OTU. If the OTU has not been classified to a given
#' taxonomic level, the entry for it is "none".
#'
#' @docType data
#' @usage data(david_stool_lineages)
#' @format a matrix with 5431 rows representing OTUs and 10 columns representing taxonomic information
#' @keywords datasets
#' @references David et al. (2014) Host lifestyle affects human microbiota on daily timescales Genome Biology vol. 15 (7):R89
#' \href{http://www.genomebiology.com/2014/15/7/R89}{Genome Biology}
#' @examples
#' # list the families with most OTUs in the David stool data
#' rev(sort(table(david_stool_lineages[,6])))[1:10]
#' # get the genus of the top abundant OTU in the stool A data set
#' data(david_stoolA_otus)
#' sorted=sort(apply(david_stoolA_otus,1,sum),index.return=TRUE, decreasing=TRUE)
#' index=which(david_stool_lineages[,1]==rownames(david_stoolA_otus)[sorted$ix[1]])
#' david_stool_lineages[index,7]
"david_stool_lineages"
