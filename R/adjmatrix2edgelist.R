#' Convert network matrix to edge list
#'
#' Convert the network matrix to edge list, in which the first column contains the regulators and second column presents the target genes. The edge weights are in the second column.
#'
#' @param adjmatrix a network matrix
#' @param cutoff threshod to cut the edge weights. Only the edges that pass the threshold will be shown in the output results.
#' @param directed logical, if \code{FLASE} the adjmatrix is transformed to symmetric matrix
#' @return A data frame of edge list.
#' @export
#' @examples
#' ##load data
#' library(RLowPCor)
#' data(gnwdata)
#' data.exp<-gnwdata$size10$ts1[,-c(1:2)]
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
