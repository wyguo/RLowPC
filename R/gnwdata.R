#' DREAM4 and GeneNetWeaver simulated datasets
#'
#' The datasets are downloaded from DREAM4 challenge (\url{https://www.synapse.org/#!Synapse:syn3049712/wiki/74630}), including size
#' 10 and size 100 time-series datasets and the corresponding refernece networks.
#' Size 10 datasets include 5 topologies and 5 experiment simulations while size 100 datasets include 5 topologies and 10 experiment
#' simulations. The details of the data can be seen in paper[1,2].
#'
#' @docType data
#'
#' @usage data(gnwdata)
#'
#' @format An object of data list
#'
#' @keywords datasets
#'
#' @references
#' [1] Schaffter T, Marbach D, Floreano D: GeneNetWeaver: In silico benchmark generation and performance profiling of network
#' inference methods. Bioinformatics 2011, 27(16):2263-2270.
#'
#' [2] Marbach D, Prill RJ, Schaffter T, Mattiussi C, Floreano D, Stolovitzky G: Revealing strengths and weaknesses of methods
#' for gene network inference. Proceedings of the National Academy of Sciences 2010, 107(14):6286-6291.
#'
#' @source \url{https://www.synapse.org/#!Synapse:syn3049712/wiki/74628}
#'
#' @examples
#' data(gnwdata)
#' data.exp<-gnwdata$size10$ts1
#' ref.edge<-gnwdata$size10$net1

"gnwdata"
