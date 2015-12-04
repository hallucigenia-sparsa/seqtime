#' David et al. stool sequencing data of subject B
#'
#' The OTUs sequenced from stool samples of subject B over time.
#' Salmonella infection (food poisoning) peaked at day 159.
#'
#' @docType data
#' @usage data(david_stoolB_otus)
#' @format a matrix with 5431 rows representing OTUs and 190 columns representing samples
#' @keywords datasets
#' @references David et al. (2014) Host lifestyle affects human microbiota on daily timescales Genome Biology vol. 15 (7):R89
#' \href{http://www.genomebiology.com/2014/15/7/R89}{Genome Biology}
#' @examples
#' data(david_stoolB_otus)
#' sorted=sort(apply(david_stoolB_otus,1,sum),decreasing=TRUE,index.return=TRUE) # sort OTUs by abundance
#' tsplot(david_stoolB_otus[sorted$ix[1:10],],type="l", header="Stool B") # display the top 10 most abundant OTUs
"david_stoolB_otus"
