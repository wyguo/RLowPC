#' Consensus network integrated from average edge ranks
#'
#' Consensus network is built from taking the average ranks of the edges in multiple network predictions.
#'
#' @param adjmatrix.list a list of inferred network matrices with same row and column ordering.
#' @param directed logical. If TRUE, the network is considered as directed. If FALSE, the upper triangular part of the symmetric matrix is used to calculate the rank matrix
#' @return \code{average.consensus} returns a network with rank average edge ranks. The weights are rescaled to 0-1 and hihger values indicate higher ranks.
#' @export
#' @examples
#' ##create two random networks
#' library(RLowPC)
#' set.seed(4)
#' net1<-abs(matrix(rnorm(16),4,4))
#' net1<-pmax(net1,t(net1))
#' diag(net1)<-0
#' set.seed(5)
#' net2<-abs(matrix(rnorm(16),4,4))
#' net2<-pmax(net2,t(net2))
#' diag(net2)<-0
#' dimnames(net1)<-dimnames(net2)<-list(letters[1:4],letters[1:4])
##put networks into list
#' net.list<-list(net1=net1,net2=net2)
##build consensus network
#' inf.consensus<-average.consensus(adjmatrix.list = net.list,directed = F)
#' adj2rankadj(net1)
#' adj2rankadj(net2)
#' inf.consensus

average.consensus<-function(adjmatrix.list,directed=F){
#   adjrank.array=abind::abind(lapply(adjmatrix.list,function(x)
#   adj2rankadj(x,directed = directed)),along=3)
#   average.rank<-apply(adjrank.array,c(1,2),mean)
  adjrank.list<-lapply(adjmatrix.list,function(x) adj2rankadj(x,directed = directed))
  average.rank<-Reduce("+", adjrank.list) / length(adjrank.list)
  average.rank<-average.rank-min(average.rank[upper.tri(average.rank)])
  average.rank<-average.rank/max(average.rank)
  diag(average.rank)<-0
  return(average.rank)
}
