#' Convert network matrix to edge list
#'
#' Convert the network matrix to edge list, in which the first column contains the regulators, the second column presents the target genes and edge weights are in the third column.
#'
#' @param adjmatrix a network matrix
#' @param cutoff threshod to cut the edge weights. Only the edges that pass the threshold will be shown in the output results.
#' @param directed logical, if \code{FLASE} the adjmatrix is transformed to symmetric matrix and the upper triangular part of the matrix is used to generate the edge list.
#' @param order logical, to order the edgelist according the edge weightes or not.
#' @return \code{adjmatrix2edgelist} returns a data frame of edge list.
#' @export
#' @examples
#' ##load data
#' library(RLowPC)
#' data(gnwdata)
#' data.exp<-gnwdata$size100$ts1[,-c(1:3)]
#' genes<-colnames(data.exp)
#' ##build correlation network
#' inf.cor<-abs(cor(data.exp))
#' diag(inf.cor)<-0
#' ##convert matrix to edge list
#' adjmatrix2edgelist(inf.cor)

adjmatrix2edgelist<-function(adjmatrix,cutoff=0,directed=F,order=T){
  if(!directed){
    adjmatrix[lower.tri(adjmatrix)]<-0
  }
  edgelist<-reshape2::melt(adjmatrix)
  colnames(edgelist)<-c('from','to','weight')
  edgelist<-edgelist[edgelist$weight>cutoff,]
  if(order)
    edgelist<-edgelist[order(edgelist$weight,decreasing = T),]
  rownames(edgelist)<-NULL
  return(edgelist)
}
