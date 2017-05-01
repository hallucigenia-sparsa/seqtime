#' @title David Stool Sequencing Data for Subject B
#'
#' @description The OTUs sequenced from stool samples of subject B over time.
#' @details Salmonella infection (food poisoning) peaked at day 159.
#'
#' @return #' @return OTUxSample matrix of size 5431 x 190
#' @keywords data, datasets
#' @docType data
#' @usage data(david_stoolB_otus)
#' @references David et al. (2014) Host lifestyle affects human microbiota on daily timescales Genome Biology vol. 15 (7):R89
#' \href{http://www.genomebiology.com/2014/15/7/R89}{Genome Biology}
#' @examples
#' data(david_stoolB_otus)
#' # Sort OTUs by abundance
#' sorted <- sort(apply(david_stoolB_otus,1,sum),
#'                  decreasing=TRUE,index.return=TRUE)
#' # Display the top 10 most abundant OTUs
#' tsplot(david_stoolB_otus[sorted$ix[1:10],], type = "l", header = "Stool B")
"david_stoolB_otus"
