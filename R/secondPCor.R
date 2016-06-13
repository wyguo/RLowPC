
#' Second order partial correlation
#'
#' Second order PC is calcuated according to the edge connections in the input adjacency matrix. The adjacency matrix can be a fully connected network (e.g. correlation matrix)
#' or a sparse network (e.g. correlation matrix after applying a edge cut-off). The correlation and p-value of a pair of nodes are measurd as (1) if they have at least 2
#' shared neighbours, a set of second order PC values and their p-values are calculated by removing all the paired combinations of neighbours. Edge correlation is deterimed by the
#' \code{quantile(p.values, probs=p.quant.prob)} p-value and the corresponding PC value. Setting \code{probs=1} indicates the maximum p-value will be selected as the significance metric.
#' (2) if they only have one shared neighbour, first order PC and the corresponding p-value is calculated as in \code{\link{firstPCor}} and (3) if they do not connect to the same nodes, the connection and p-value are calculated from the zero order PC (the general correlation).
#' @param adjmatrix an adjacency matrix with 0 and 1 numerical values to indicate edge connections.
#' @param data.exp gene expression data matrix with variables in columns and samples in rows.
#' @param method estimators of correlation. Options are "pearson" and "spearman"
#' @param p.quant.prob distribution probability of p-value quantile cut-off. It is used to determine which p-value of all the range of p-values for a pair of nodes is selected to determine
#'  the edge correlation.
#'
#' @return A list includes 3 \code{dataframe} of results. \code{cor} is the results of first order PC assigned to the input adjacency matrix. \code{pval} is the
#' p-value results and \code{adj} is the input adjacency matrix.
#'
#' @export




secondPCor<-function(adjmatrix,data.exp,method='pearson',p.quant.prob=1,progressbar=T){
  adjmatrix[adjmatrix>0]<-1
  ##calculate first order partial correlation
  first.stat<-firstPCor(adjmatrix = adjmatrix,data.exp = data.exp,method = method,p.quant.prob = p.quant.prob,progressbar = progressbar)
  first.cor<-first.stat$cor
  first.p<-first.stat$pval
  ##extract edge list from adjacency matrix
  adjmatrix[adjmatrix>0]<-1
  diag(adjmatrix)<-0
  edgelist<-adjmatrix2edgelist(adjmatrix,directed = T)
  ##calculate second order PC
  second.cor<-second.p<-first.cor*0
  message('Calulate second order partial correlation...')
  if(progressbar)
    pb <- txtProgressBar(min = 0, max =nrow(edgelist), style = 3)
  for(i in 1:nrow(edgelist)){
    Sys.sleep(0)
    node1<-as.vector(edgelist$from[i])
    node2<-as.vector(edgelist$to[i])
    sub.edge<-shared.neighbour(node1,node2,edgelist,verbose = F)
    if(dim(sub.edge)[1]<4){
      second.cor[node1,node2]<-first.cor[node1,node2]
      second.p[node1,node2]<-first.p[node1,node2]
    } else {
      regulators<-as.vector(unique(sub.edge$from))
      targets<-c(node1,node2)
      reg.group<-combn(x = regulators,2)
      cor.v<-p.v<-vector()
      for(j in 1:ncol(reg.group)){
        sub.reg<-reg.group[,j]
        second.stat<-ppcor::pcor(data.exp[,c(sub.reg,targets)],method = method)
        cor.v<-c(cor.v,second.stat$estimate[node1,node2])
        p.v<-c(p.v,second.stat$p.value[node1,node2])
      }
      second.cor[node1,node2]<-cor.v[p.v==quantile(p.v,probs = p.quant.prob)]
      second.p[node1,node2]<-quantile(p.v,probs = p.quant.prob)
    }
    if(progressbar)
      setTxtProgressBar(pb, i)
  }
  if(progressbar)
    close(pb)
  results<-list(cor=second.cor,pval=second.p,adj=adjmatrix)
  return(results)
}

