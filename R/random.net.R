#' Resort datasets and generate random networks
#'
#' The function is used to generate random datasets or random networks
#' @details
#' (a) Input a gene expression data frame with genes in columns and samples in rows, the function resorts the sample points for each genes in each experiments.
#' (b) Input a string vector of gene names, the function generates random network with random edge weightes.
#'
#' @param input a gene expression data frame or a string vector of gene names.
#' @param nexp a numeric number to indicate how many experiments in the datasets. The gene expression datasets are resorted in both experiment-wise and gene-wise.
#' @param type a character string to indicate the resort datasets ("data") or to generate random networks ("network").
#' @param directed logical, if \code{FALSE}, the random network generated from gene names will be converted to symmetric matrix.
#' @return \code{random.net} returns a randomly resorted expression data frame or a random network matrix.
#' @export
random.net<-function(input,nexp=NA,type='data',directed=F){
  if(type=='data'){
    if(!is.data.frame(input))
      stop('The input does not match to the "type" to generate random networks.')
    if(is.na(nexp))
      stop('Please provide how many experiments in the time-series datasets.')
    data4random<-data.frame(id=rep(1:nexp,each=nrow(input)/nexp),input)
    result<-plyr::ddply(data4random,.(id),colwise(sample))[,-1]
  }
  if(type=='network'){
    if(!is.vector(input))
      stop('The input does not match to the "type" to generate random networks.')
    inf.random<-matrix(runif(length(input)^2),length(input),length(input))
    dimnames(inf.random)<-list(input,input)
    diag(inf.random)<-0
    if(!directed)
      inf.random<-pmax(inf.random,t(inf.random))
    result<-inf.random
  }
  return(result)
}
