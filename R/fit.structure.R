#' Calculate structural fitness of sub-networks
#'
#' The function is to used to calcuate the scale-free or modularity fitness of top weighted sub-networks.The scale-free fitness is calcuated as
#' \eqn{-R^2\times slope} of the fit index of the model. Please refer to the package \code{\link[WGCNA]{WGCNA}} [1,2] for details.The fitness of modularity is computed
#' using \code{\link[igraph]{igraph}} package [3].

#' @param net a network matrix or data frame of edge list.
#' @param top.vector a numeric vector of top weighted edge number, e.g. \code{top.vector=seq(100,500,by=10)}.
#' @param method a character string used to calculate the structural fitness. Options are "scale-free" and "modularity"
#' @return a vector of fitness
#' @references
#' [1] Langfelder P, Horvath S: WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559.
#'
#' [2] Zhao W, Langfelder P, Fuller T, Dong J, Li A, Hovarth S: Weighted gene coexpression network analysis: state of the art. J Biopharm Stat 2010, 20(2):281-300.
#'
#' [3] Csardi G, Nepusz T: The igraph software package for complex network research. InterJournal 2006, Complex Systems:1695.
#'
#' @export
#' @examples
#'library(networkBMA)
#'library(RLowPCor)
#'library(minet)
#'data(dream4)
#'ref.edge<-dream4gold100[[1]]
#'top.vector<-seq(100,300,by=10)
#'top.network(ref.edge,top.vector=top.vector,method='scale-free')

fit.structure<-function(net,top.vector,method='scale-free'){
  if(ncol(net)!=nrow(net) & ncol(net)!=3)
    stop('Please provide a network matrix or edge list')
  if(ncol(net)==nrow(net))
    net<-adjmatrix2edgelist(net,directed = T)
  net<-net[order(net[,3],decreasing = T),]
  top.fitness<-vector()
  if('scale-free' %in% method){
    for(top in top.vector){
      sub.net<-net[1:top,]
      sub.degree<-node.degree(sub.net,directed = T)
      fitness<-WGCNA::scaleFreeFitIndex(sub.degree)
      fitness<-as.numeric(fitness[1])
      top.fitness<-c(top.fitness,fitness)
    }
  }
  if('modularity' %in% method){
    for(top in top.vector){
      sub.net<-net[1:top,]
      g<-igraph::graph.data.frame(sub.net,directed = T)
      wtc <- igraph::walktrap.community(g)
      top.fitness<-c(top.fitness, igraph::modularity(wtc))
    }
  }
  names(top.fitness)<-top.vector
  return(top.fitness)
}
