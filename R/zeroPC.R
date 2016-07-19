
#' Zero order partial correlation
#'
#'The correlation and p-value of correlation.
#' @param data.exp gene expression data.
#' @param method "pearson", "spearman" or "kendall" method is used to calcuate the correlation values.
#'
#' @examples
#' ##load data
#' data(gnwdata)
#' ##ts expression
#' data.exp<-gnwdata$size100$ts1[,-c(1:3)]
#' inf.zeroPC<-zeroPC(data.exp)
#' head(inf.zeroPC)
#' ##refernece network
#'
#' @return \code{zeroPC} returns a data frame with columns of regulators, target genes, edge correlation weigthes,
#' p-values, FRD and connection probability.
#'
#' @export

zeroPC<-function(data.exp,method='pearson'){
  message('Calulate zero order partial correlation (the general correlation)...')
  inf.cor<-abs(cor(data.exp,method = method))
  diag(inf.cor)<-0
  inf.edge<-adjmatrix2edgelist(inf.cor,directed = F,order = T)
  ##calculate statistics, pval, fdr and prob
  fdr.frame<-cor2statistics(cor.vector = abs(inf.edge$weight),plot=F,verbose=F)
  inf.edge$pval = fdr.frame$pval
  inf.edge$qval = fdr.frame$qval
  inf.edge$prob = 1 - fdr.frame$lfdr
  rownames(inf.edge)<-NULL
  return(inf.edge)
}
