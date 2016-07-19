#' Convert edge list to network matrix
#'
#' The function is to convert edge list to network matrix.
#'
#' @param edgelist a data frame of edge list of a network, in which columns are regulators, targets and edge weights.
#' @param genes gene names to name the rows and columns of the output network matrix.
#' @param directed logical, to create directed or undirected (symmetric) network matrix.
#' @param cutoff the threshold to cut the edge list.
#' @return \code{edgelist2adjmatrix} returns a network matrix.
#' @export
#' @examples
#' ##load data
#' library(RLowPC)
#' data(gnwdata)
#' data.exp<-gnwdata$size100$ts1[,-c(1:3)]
#' genes<-colnames(data.exp)
#' ref.edge<-gnwdata$size100$net1
#' ref.edge[,3]<-1
#' ref.adj<-edgelist2adjmatrix(ref.edge,genes,directed=F)

edgelist2adjmatrix<-function(edgelist,genes,cutoff=0,directed=F){
  colnames(edgelist)<-c('from','to','weight')
  edgelist<-edgelist[edgelist$weight>cutoff,]
  adjmatrix<-matrix(0,ncol=length(genes),nrow=length(genes))
  dimnames(adjmatrix)<-list(genes,genes)
  for(i in 1:nrow(edgelist)){
    adjmatrix[as.vector(edgelist[i,1]),as.vector(edgelist[i,2])]<-edgelist[i,3]
  }
  if(!directed){
    adjmatrix<-pmax(adjmatrix,t(adjmatrix))
  }
  return(adjmatrix)

}
