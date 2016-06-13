#' Convert correlation matrix to MI matrix
#'
#' The correlation \eqn{\rho_{ij}} between node \eqn{i} and \eqn{j} is converted into MI by using equation [1]:
#' \deqn{MI=-\frac{1}{2} (1-\rho_{ij}^2)}.
#'
#' @param cor.adj a correlation matrix
#' @references
#' [1] Meyer PE, Kontos K, Lafitte F, Bontempi G: Information-theoretic inference of large transcriptional
#' regulatory networks. EURASIP J Bioinform Syst Biol 2007:79879.
#' @return A mutual information matrix
#' @export
cor2mi=function(cor.adj){
  mi=-0.5*log(1-x^2)
  return(mi)
}
