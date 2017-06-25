#' @title Caporaso Stool Sequencing Data for Subject M3 on level 6
#'
#' @description Mostly genus-level counts (L6) from stool samples of subject M3 over time.
#' @details Subject M3 was followed for 15 months. V4 16S region was sequenced. A counter was attached to genera appearing more than once in the count table.
#' @return TaxonxSample matrix of size 335 x 332
#' @docType data
#' @usage data(caporaso_M3FecesL6)
#' @keywords data, datasets
#' @references
#' Caporaso et al. (2011) Moving pictures of the human microbiome Genome Biology vol. 12:R50
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r50}{Genome Biology}
#' @examples
#' data(caporaso_M3FecesL6)
#' # sort taxa by abundance and print top 10 most abundant taxa
#' sorted <- sort(apply(caporaso_M3FecesL6,1,sum),decreasing = TRUE, index.return = TRUE)
#' rownames(caporaso_M3FecesL6)[sorted$ix[1:10]]
"caporaso_M3FecesL6"
