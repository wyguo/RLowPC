#' Extract sub-edge list with shared neihbours
#'
#' The function is used to extract the shared neighbours that both connect to a pair of candidate nodes.
#'
#' @param node1 a character string of the name of candidate node 1
#' @param node2 a character string of the name of candidate node 2
#' @param edgelist the full edge list used to extract the sub-edge list. First column are regulators, second column are target genes and the last column are weigths of the connections.
#' @return Sub-edge list
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
