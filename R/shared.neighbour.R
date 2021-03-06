#' Extract shared neihbour genes
#'
#' The function is used to extract the shared neighbour genes that both connect to a pair of candidate nodes.
#'
#' @param node1 a character string of the name of candidate node 1.
#' @param node2 a character string of the name of candidate node 2.
#' @param edgelist the edge list of a network from which the neighbour genes are extracted for node1 and node2.
#' @return \code{shared.neighbour} returns a sub edgelist.
#' @export
shared.neighbour<-function(node1, node2, edgelist,verbose=T){
  colnames(edgelist)<-c('from','to','weight')
  sub.edgelist<-edgelist[which(edgelist$to %in% node1 | edgelist$to %in% node2),]
  ##filter the neigbours not shared by node1 and node2
  filter.gene<-table(sub.edgelist[,1])
  filter.gene<-names(filter.gene)[filter.gene==1]
  if(length(filter.gene)!=0){
    sub.edgelist<-sub.edgelist[-which(sub.edgelist[,1] %in% filter.gene),]
  }
  if(nrow(sub.edgelist)==0 & verbose)
    message(paste0("There is no shared neighbour for ",node1," and ", node2))
  return(sub.edgelist)
}
