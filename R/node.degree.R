#' Calculate node degree of gene network
#'
#' @param net input network matrix or edge list with three columns, i.e. regulators, targets and weights.
#' @param directed logical, if \code{FALSE}, the input network is treated as undirected.
#' @return a vector of node degrees
#' @export
node.degree<-function(net,directed=F){
  if(ncol(net)!=nrow(net) & ncol(net)!=3)
    stop('Please provide a network matrix or edge list')
  if(ncol(net)==3){
    net<-net[net[,3]>0,]
    degree<-table(unlist(net[,1:2]))
  }
  if(ncol(net)==nrow(net)){
  net[net>0]<-1
  degree<-colSums(net)+rowSums(net)
  }
  if(!directed)
    degree<-degree/2
  return(degree)
}
