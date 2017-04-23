#' @title Constructor for S3 noisetypes class
#'
#' @description Construct a noisetypes object.
#'
#' @param white indices of matrix rows with white noise
#' @param pink indices of matrix rows with pink noise
#' @param brown indices of matrix rows with brown noise
#' @param black indices of matrix rows with black noise
#' @return noisetypes class
#' @export

noisetypes<-function(white=0, pink=0, brown=0, black=0){
  noisetypes=list(white,pink,brown,black)
  names(noisetypes)=c("white","pink","brown","black")
  attr(noisetypes, "class") <- "noisetypes"
  return(noisetypes)
}
