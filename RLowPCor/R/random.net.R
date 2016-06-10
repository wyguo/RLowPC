#' Generate random network
#'
#' Random networks are generated from node size or from resort odering of samples in gene expression dataset.
#' @param input a gene expression \code{dataframe} with variables in columns and samples in rows or a vector of node names.
#' @param nblock a numeric number to indicate how many blocks that consist the datasets. For example, to tell how many replicates or experiments are in the datasets.
#' The datasets to generate random network will be resorted in block-wise.
#' @param from a character string to tell the random network generation is based on nodes size or dataset.Options are "node" and "data".
#' @param directed logical, if \code{FALSE}, the random network generated from nodes will be converted to symmetric matrix.
#' @return The resorted datasets if \code{from='data'} or a random network matrix, the edge weights of which are in \eqn{[0,1]} if \code{from='node'}.
#' @export
random.net<-function(input,nblock,from='data',directed=F){
  if(from=='data'){
    if(!is.data.frame(input))
      stop('The input does not match to the way to generate random networks.')
    data4random<-data.frame(id=rep(1:nblock,each=nrow(input)/nblock),input)
    result<-plyr::ddply(data4random,.(id),colwise(sample))[,-1]
  }
  if(from=='node'){
    if(!is.vector(input))
      stop('The input does not match to the way to generate random networks.')
    inf.random<-matrix(runif(length(input)^2),length(input),length(input))
    dimnames(inf.random)<-list(input,input)
    diag(inf.random)<-0
    if(!directed)
      inf.random<-pmax(inf.random,t(inf.random))
    result<-inf.random
  }
  return(result)
}
