#' Convert network matrix to edge list
#'
#' The inferred network matrix is converted to edge list
#'
#' @param adjmatrix a network matrix
#' @param cutoff threshod to cut the edge list
#' @param directed logical, if \code{FLASE} the adjmatrix is transformed to symmetric matrix
#' @return a edge list
#' @export
#' @examples
#' ##load data
#' library(networkBMA)
#' library(RLowPCor)
#' data(dream4)
#' data.exp<-dream4ts10[[1]][,-c(1:2)]
#' genes<-colnames(data.exp)
#' ##build correlation network
#' inf.cor<-abs(cor(data.exp))
#' diag(inf.cor)<-0
#' ##convert matrix to edge list
#' adjmatrix2edgelist(inf.cor)

adjmatrix2edgelist<-function(adjmatrix,cutoff=0,directed=F){
  if(!directed){
    adjmatrix<-pmax(adjmatrix,t(adjmatrix))
  }
  edgelist<-reshape2::melt(adjmatrix)
  colnames(edgelist)<-c('from','to','weight')
  edgelist<-edgelist[order(edgelist$weight,decreasing = T),]
  edgelist<-edgelist[edgelist$weight>cutoff,]
  rownames(edgelist)<-NULL
  return(edgelist)
}
