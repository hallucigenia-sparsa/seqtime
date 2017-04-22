#' @title Constructor for S3 noisetypes class
#'
#' @description Construct a noisetypes object.
#'
#' @param white percentage of white noise rows in matrix
#' @param pink percentage of pink noise rows in matrix
#' @param brown percentage of brown noise rows in matrix
#' @param black percentage of black noise rows in matrix
#' @return noisetypes class
#' @export

noisetypes<-function(white=0, pink=0, brown=0, black=0){
  nonclass=100-(white+pink+brown+black)
  noisetypes=list(white,pink,brown,black,nonclass)
  names(noisetypes)=c("white","pink","brown","black","nonclass")
  attr(noisetypes, "class") <- "noisetypes"
  return(noisetypes)
}
