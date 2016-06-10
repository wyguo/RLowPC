#' Remove the zero rows and columns of data frame
#'
#' The function is used to remove the zero rows and zero columns in network matrix
#' @param adjmatrix a network matrix
#' @return A new network matrix
#' @export
remove.zero<-function(adjmatrix){
  adjmatrix<-adjmatrix[rowSums(adjmatrix)>0,]
  adjmatrix<-adjmatrix[,colSums(adjmatrix)>0]
  return(adjmatrix)
}
