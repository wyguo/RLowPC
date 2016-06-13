#' Convert network matrix to network rank matrix
#'
#' @param adjmatrix A network matrix
#' @param directed Logical. If TRUE, the network is considered as directed. If FALSE, the upper triangular part of the symmetric network matrix is used to calcuate the rank matrix.
#' @return A network matrix with rank weightes.
#' @export
adj2rankadj<-function(adjmatrix,directed=F){
  adjmatrix.rank<-adjmatrix-adjmatrix
  if(!directed){
    values<-rank(adjmatrix[upper.tri(adjmatrix)])
    adjmatrix.rank[upper.tri(adjmatrix.rank)]<-values
    adjmatrix.rank<-pmax(adjmatrix.rank,t(adjmatrix.rank))
  }
  if(directed){
    diag(adjmatrix)<-NA
    values<-rank(adjmatrix)
    adjmatrix.rank[]=values
    diag(adjmatrix.rank)<-0
  }
  return(adjmatrix.rank)
}
