#' Convert edge list to network matrix
#'
#' The function is to convert edge list to network matrix.
#'
#' @param edgelist a data frame of edge list of a network, in which columns are regulators, targets and edge weights.
#' @param genes gene names to name the rows and columns of the output network matrix.
#' @param directed logical, to create directed or undirected () network matrix.
#' @param cutoff the threshold to cut the edge list
#' @return a network matrix
#' @export
#' @examples
#'  library(networkBMA)
#'  library(RLowPCor)
#'  ##load DREAM4 size100_1 datasets
#'  data(dream4)
#'  data.exp<-dream4ts10[[1]][,-c(1:2)]
#'  genes<-colnames(data.exp)
#'  ref.edge<-dream4gold10[[1]]
#'  ref.adj<-edgelist2adjmatrix(ref.edge,genes)

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
