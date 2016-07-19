#' Evaluation of inferred gene networks
#'
#' The inferred network is evaluated by comparing to the reference network. The output is a tables of TP, FP, TN and FP with different edge weight cut-offs [1].
#'
#' @param inf.adj the inferred network matrix. Column names and row names match to the reference network.
#' @param ref.adj the reference network matrix with 1 inidating connected edge and 0 unconnected edge.
#' @param directed logical, to compare as directed or undirected networks. In a undirected network, only the upper triangular of the network matrix is used for evaluation.
#' @references
#' [1] Meyer PE, Lafitte F, Bontempi G: minet: A R/Bioconductor package for inferring large transcriptional networks using mutual information.
#' BMC Bioinformatics 2008, 9:461.
#' @export
#' @examples
#' ##load library
#' library(RLowPC)
#' ##load data
#' data(gnwdata)
#' data.exp<-gnwdata$size100$ts1[,-c(1:3)]
#' genes<-colnames(data.exp)
#' ref.edge<-gnwdata$size100$net1
#' ref.edge[,3]<-1
#' ref.adj<-edgelist2adjmatrix(ref.edge,genes)
#' inf.cor<-abs(cor(data.exp))
#' diag(inf.cor)<-0
#' table.cor<-table.evaluate(inf.adj = inf.cor,ref.adj = ref.adj)
#' head(table.cor)

table.evaluate<-function(inf.adj,ref.adj,directed=F){
  if(directed){
    table.output<-minet::validate(inf.adj,ref.adj)
    table.zero<-table.output[table.output$thrsh==0,]
    table.output<-rbind(table.output[table.output$thrsh>0,],table.zero[1:10,])
  } else {
    inf.adj<-pmax(inf.adj,t(inf.adj))
    ref.adj<-pmax(ref.adj,t(ref.adj))
    ref.adj[lower.tri(ref.adj)]<-inf.adj[lower.tri(inf.adj)]<-0
    table.output<-minet::validate(inf.adj,ref.adj)
    tn.remove<-(ncol(ref.adj)^2-ncol(ref.adj))/2+ncol(ref.adj)
    table.output$tn<-table.output$tn-tn.remove
    table.output<-table.output[table.output$tn>=0,]
    table.zero<-table.output[table.output$thrsh==0,]
    table.output<-rbind(table.output[table.output$thrsh>0,],table.zero[1:10,])
    rownames(table.output)<-NULL
  }
  table.output<-na.omit(table.output)
  return(table.output)
}

