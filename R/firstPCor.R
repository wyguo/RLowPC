
#' First order partial correlation
#'
#' First order PC is calcuated according to the edge connections in the input adjacency matrix. The adjacency matrix can be a fully connected network (e.g. correlation matrix)
#' or a sparse network (e.g. correlation matrix after applying a edge cut-off). The correlation and p-value of a pair of nodes are measurd as (1) if they have
#' shared neighbours, a set of first order PC values and their p-values are calculated by removing the neighbours one by one. Edge correlation is deterimed by the
#' \code{quantile(p.values, probs=p.quant.prob)} p-value and the corresponding PC value. Setting \code{probs=1} indicates the maximum p-value will be selected as the significance metric.
#' (2) if they do not connect to the same nodes, the connection and p-value are calculated from the zero order PC (the general correlation).
#' @param adjmatrix an adjacency matrix with 0 and 1 numerical values to indicate edge connections.
#' @param data.exp gene expression data matrix with variables in columns and samples in rows.
#' @param method estimators of correlation. Options are "pearson" and "spearman"
#' @param p.quant.prob distribution probability of p-value quantile cut-off. It is used to determine which p-value of all the range of p-values for a pair of nodes is selected to determine
#'  the edge correlation.
#' @param progressbar logical. If TRUE, a progressbar will show to indicate the code runing percentage.
#' @examples
#'  library(RLowPCor)
#'  ## load DREAM4 data
#'  data(gnwdata)
#'  data.exp<-gnwdata$size100$ts1[,-c(1:2)]
#'  genes<-colnames(data.exp)
#'  ref.edge<-gnwdata$size100$net1
#'  ref.adj<-edgelist2adjmatrix(ref.edge,genes)
#'  ## infer correlation network
#'  inf.cor<-abs(cor(data.exp))
#'  diag(inf.cor)<-0
#'  ##  infer shrink partial correlation network
#'  inf.pcor<-abs(pcor.shrink(data.exp))[1:ncol(data.exp),1:ncol(data.exp)]
#'  diag(inf.pcor)<-0
#'  adjmatrix<-inf.cor
#'  adjmatrix[adjmatrix<0.2]<-0
#'  ## infer first order partial correlation network
#'  inf.firstpcor<-abs(firstPCor(adjmatrix = adjmatrix,data.exp)$cor)
#'  ##  infer second order partial correlation network
#'  inf.firstpcor<-abs(secondPCor(adjmatrix = adjmatrix,data.exp)$cor)
#'
#'  table.firstpcor<-table.evaluate(inf.firstpcor,ref.adj)
#'  table.secondpcor<-table.evaluate(inf.secondpcor,ref.adj)
#'  table.cor<-table.evaluate(inf.cor,ref.adj)
#'  table.pcor<-table.evaluate(inf.pcor,ref.adj)
#'  X11()
#'  plotAUC(list(cor=table.cor,firstpcor=table.firstpcor,secondpcor=table.secondpcor,shrink.pcor=table.pcor),lwd=2)

#'
#' @return A list includes 3 \code{dataframe} of results. \code{cor} is the results of first order PC assigned to the input adjacency matrix. \code{pval} is the
#' p-value results and \code{adj} is the input adjacency matrix.
#'
#' @export

firstPCor<-function(adjmatrix,data.exp,method='pearson',p.quant.prob=1,progressbar=T){
  adjmatrix[adjmatrix>0]<-1
  diag(adjmatrix)<-0
  edgelist<-adjmatrix2edgelist(adjmatrix,directed = T)

  ###calculate zero order correlation
  zero.stat=Hmisc::rcorr(as.matrix(data.exp),type = method)
  zero.cor<-zero.stat$r
  diag(zero.cor)<-0
  zero.p<-zero.stat$P
  diag(zero.p)<-0
  ####
  first.cor<-first.p<-adjmatrix*0
  message('Calulate first order partial correlation...')
  if(progressbar)
    pb <- txtProgressBar(min = 0, max =nrow(edgelist), style = 3)
  for(i in 1:nrow(edgelist)){
    Sys.sleep(0)
    node1<-as.vector(edgelist$from[i])
    node2<-as.vector(edgelist$to[i])
    sub.edge<-shared.neighbour(node1,node2,edgelist,verbose = F)
    if(dim(sub.edge)[1]==0){
      first.cor[node1,node2]<-zero.cor[node1,node2]
      first.p[node1,node2]<-zero.p[node1,node2]
    } else {
      regulators<-as.vector(unique(sub.edge$from))
      targets<-c(node1,node2)
      cor.v<-p.v<-vector()
      for(sub.reg in regulators){
        first.stat<-ppcor::pcor(data.exp[,c(sub.reg,targets)],method = method)
        cor.v<-c(cor.v,first.stat$estimate[node1,node2])
        p.v<-c(p.v,first.stat$p.value[node1,node2])
      }
      first.cor[node1,node2]<-cor.v[p.v==quantile(p.v,probs = p.quant.prob)]
      first.p[node1,node2]<-quantile(p.v,probs = p.quant.prob)
    }
    if(progressbar)
    setTxtProgressBar(pb, i)
  }
  if(progressbar)
    close(pb)
  results<-list(cor=first.cor,pval=first.p,adj=adjmatrix)
  return(results)
}

