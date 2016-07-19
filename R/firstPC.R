#' First order partial correlation
#'
#' First order PC is calcuated according to the edge connections in the input edgelist. The correlation of a pair of genes are calcuated as (a) correlation if they do not connect
#' to the same neighbour genes and (b) the maximum of first order PC values by removing the controls one by one if they connect to shared neighbour genes. After that, the
#' \code{\link{cor2statistics}} is used to calculate the p-value, FDR and connection probability from the output correlation.
#' @param edgelist an edge list for fisrt order PC calculation. The PC values are calcuated over each pair of genes in the gene list.
#' @param controlist list of vectors of neighbour genes for each pair of genes in \code{edgelist}. The vector order in the list must match to the order of corresponding paired genes in the
#' input \code{edgelist}. For example the i element in \code{controlist} contains the shared neighbours of gene \code{edgelist[i,1]} and gene \code{edgelist[i,2]}.
#' @param data.exp gene expression data matrix with genes in columns and samples in rows.
#' @param method estimators of correlation. Options are "pearson", "spearman" and "kendall".
#' @param logical. If TRUE, a progressbar will show to indicate the code runing percentage.
#'
#' @return \code{firstPC} returns a data frame with columns of regulators, target genes, edge correlation weigthes, p-values, FRD and connection probability.
#' @seealso \code{\link{cor2statistics}}, \code{\link{RLowPC}}
#'
#' @export

firstPC<-function(data.exp,edgelist,controlist=NULL,method='pearson',progressbar=T){
  edgelist.new<-edgelist
  ####
  t1 <- Sys.time()
  ##zero ordr PC
  inf.zeroPC<-zeroPC(data.exp=data.exp,method = method)

  message(paste0('Calulate first order partial correlation for ',nrow(edgelist),' pairs of genes ...'))
  if(progressbar)
    pb <- txtProgressBar(min = 0, max =nrow(edgelist), style = 3)
  for(i in 1:nrow(edgelist)){
    Sys.sleep(0)
    node1<-as.vector(edgelist$from[i])
    node2<-as.vector(edgelist$to[i])
    sub.edge<-shared.neighbour(node1,node2,edgelist,verbose = F)
    if(dim(sub.edge)[1]==0){
      edgelist.new[i,3]<-inf.zeroPC[inf.zeroPC$from==node1 &inf.zeroPC$to==node2,]$weight
    } else {
      if(is.null(controlist))
        regulators<-as.vector(unique(sub.edge$from))
      if(!is.null(controlist))
        regulators<-controlist[[i]]
      if(node1 %in% regulators)
        regulators<-regulators[-which(regulators %in% node1)]
      if(node2 %in% regulators)
        regulators<-regulators[-which(regulators %in% node2)]
      targets<-c(node1,node2)
      cor.v<-p.v<-vector()
      for(sub.reg in regulators){
        first.stat<-ppcor::pcor(data.exp[,c(sub.reg,targets)],method = method)
        first.stat.cor<-first.stat$estimate
        dimnames(first.stat.cor)<-list(c(sub.reg,targets),c(sub.reg,targets))
        cor.v<-c(cor.v,abs(first.stat.cor[node1,node2]))
      }

      edgelist.new[i,3]<-cor.v[abs(cor.v)==max(abs(cor.v))][1]
    }
    if(progressbar)
    setTxtProgressBar(pb, i)
  }
  if(progressbar)
    close(pb)
  edgelist.new<-edgelist.new[order(abs(edgelist.new$weight),decreasing = T),]
  fdr.frame<-cor2statistics(cor.vector = abs(edgelist.new$weight),plot=F,verbose=F)
  edgelist.new$pval = fdr.frame$pval
  edgelist.new$qval = fdr.frame$qval
  edgelist.new$prob = 1 - fdr.frame$lfdr
  rownames(edgelist.new)<-NULL
  t2 <- Sys.time()
  message(paste0('Done! Time taken:',round(as.numeric(difftime(t2,t1)),4),' ',units(difftime(t2,t1))))
  return(edgelist.new)
}

