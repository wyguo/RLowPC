#' Select sub-network with top weigthed edges
#'
#' @param net a netork matrix or edge list
#' @param top a numeric number to indicate the top weigthed edges, e.g. \code{top=1000}
#' @return A subset of network matrix or edge list
#' @export


top.network<-function(net,top=1000){
  if(isSymmetric(net))
    net[lower.tri(net)]<-0
  if(nrow(net)==ncol(net))
    edge<-adjmatrix2edgelist(net,directed = T)
  sub.net<-edge[1:top,]
  if(nrow(net)==ncol(net)){
    net.new<-edgelist2adjmatrix(sub.net,genes=colnames(net),directed = T)
  }
}
