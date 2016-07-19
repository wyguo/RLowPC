#' Convert network matrix to network rank matrix
#'
#' The edge weights in the network matrix will be converted to their ranks. A high rank value indcates a high edge connection.
#'
#' @param adjmatrix A network matrix
#' @param directed Logical. If TRUE, the network is considered as directed. If FALSE, the upper triangular part of the symmetric network matrix is used to calcuate the rank matrix.
#' @return \code{adj2rankadj} returns a network matrix with ranks of the edge weightes.
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
