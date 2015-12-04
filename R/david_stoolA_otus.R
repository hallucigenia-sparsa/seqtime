#' David et al. stool sequencing data of subject A
#'
#' The OTUs sequenced from stool samples of subject A over time.
#' Subject A traveled from day 71 to day 122 and contracted diarrhea on
#' days 80 to 85 and 104 to 113.
#'
#' @docType data
#' @usage data(david_stoolA_otus)
#' @format a matrix with 5431 rows representing OTUs and 329 columns representing samples
#' @keywords datasets
#' @references David et al. (2014) Host lifestyle affects human microbiota on daily timescales Genome Biology vol. 15 (7):R89
#' \href{http://www.genomebiology.com/2014/15/7/R89}{Genome Biology}
#' @examples
#' data(david_stoolA_otus)
#' sorted=sort(apply(david_stoolA_otus,1,sum),decreasing=TRUE,index.return=TRUE) # sort OTUs by abundance
#' tsplot(david_stoolA_otus[sorted$ix[1:10],],type="l", header="Stool A") # display the top 10 most abundant OTUs
"david_stoolA_otus"
