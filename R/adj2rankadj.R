#' Convert network matrix to network rank matrix
#'
#' @param adjmatrix A network matrix
#' @param directed Logical. If TRUE, the network is considered as directed. If FALSE, the upper triangular part of the matrix is used to calcuate the rank matrix
#' @return an network connection rank matrix
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
