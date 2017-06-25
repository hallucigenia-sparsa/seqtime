#' @title David Stool Sequencing Data for Subject A
#'
#' @description The OTUs sequenced from stool samples of subject A over time.
#' @details Subject A traveled from day 71 to day 122 and contracted diarrhea on
#' days 80 to 85 and 104 to 113.
#' @return OTUxSample matrix of size 5431 x 329.
#' @docType data
#' @usage data(david_stoolA_otus)
#' @keywords data, datasets
#' @references
#'    David et al. (2014) Host lifestyle affects human microbiota on daily timescales Genome Biology vol. 15 (7):R89
#' \href{http://www.genomebiology.com/2014/15/7/R89}{Genome Biology}
#' @examples
#' data(david_stoolA_otus)
#'   # sort OTUs by abundance
#'   sorted <- sort(apply(david_stoolA_otus,1,sum),
#'                  decreasing = TRUE, index.return = TRUE)
#' # display the top 10 most abundant OTUs
#' tsplot(david_stoolA_otus[sorted$ix[1:10],], header = "Stool A")
"david_stoolA_otus"
