#' Convert correlation to statistical significance measures.
#'
#'\code{cor2statistics} convert a correlation vector to p-values, false discovery rates (FDRs) and probability of connections.
#' @param cor.vector a correlation vector
#' @param ... parameters used in \code{\link{fdrtool}}.
#'
#' @return \code{cor2statistics} returns a list includes correlation matrix and p-value matrix.
#' @seealso \code{\link{fdrtool}}
#'
#' @export


cor2statistics<-function(cor.vector,...){
  fdr.frame = fdrtool::fdrtool(abs(cor.vector), statistic = "correlation",...)
  return(fdr.frame)
}
