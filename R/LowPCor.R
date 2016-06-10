
#' Calculate low order partial correlation
#'
#' Calculate up to second order partial correlation (PC).The zero order PC is the general correlation. The first order PC is
#' calcualted based on zero PC matrix after applying correlation significance testing. In the zero PC matrix, the less significant edges are
#' removed with a p-value cut-off. The second order PC is calcuated based on first order PC matrix after significance testing. For all the combinations of first order
#' PC or second order PC of a pair of node, the \code{quantile(p.values,probs=p.quant.prob)} quantile of p-values is selected as the
#' Such procedures reduce the computational cost, especially in large scale networks. To use fully connected network to calculate the first and second order PC,
#' users can set the p-value cut-off to \code{p.cut=1.1}.
#'
#' @param adjmatrix A network matrix
#' @param directed Logical. If TRUE, the network is considered as directed. If FALSE, the upper triangular part of the matrix is used to calcuate the rank matrix
#' @return an network connection rank matrix
#' @export




LowPCor<-function(data.exp,p.cut=0.05,up2=2,p.quant.prob=1,method='pearson',progressbar=T){
  results.list<-list()
  if(up2>=0){

    zero.stat=Hmisc::rcorr(as.matrix(data.exp),type = method)
    zero.cor<-zero.stat$r
    diag(zero.cor)<-0
    zero.p<-zero.stat$P
    diag(zero.p)<-0
    ##adjmatrix
    zero.adj<-zero.p
    zero.adj[zero.adj>=p.cut]<-0
    diag(zero.adj)<-0
    zero.adj[zero.adj>0]<-1
    ##
    results<-list(cor=zero.cor,pval=zero.p,adj=zero.adj)
    results.list<-c(results.list,setNames(list(results),'up2zero'))
  }
  if(up2>=1){
    zero2one.edge<-adjmatrix2edgelist(zero.adj)
    first.cor<-first.p<-zero.cor*0
    message('Calulate first order partial correlation...')
    if(progressbar)
      pb <- txtProgressBar(min = 0, max =nrow(zero2one.edge), style = 3)
    for(i in 1:nrow(zero2one.edge)){
      Sys.sleep(0)
      node1<-as.vector(zero2one.edge$from[i])
      node2<-as.vector(zero2one.edge$to[i])
      sub.edge<-shared.neighbour(node1,node2,zero2one.edge,verbose = F)
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
    ##adjmatrix
    first.adj<-first.p
    first.adj[first.adj>=p.cut]<-0
    diag(first.adj)<-0
    first.adj[first.adj>0]<-1
    ##
    results<-list(cor=first.cor,pval=first.p,adj=first.adj)
    results.list<-c(results.list,setNames(list(results),'up2first'))
  }
  if(up2>=2){
    first2second.edge<-adjmatrix2edgelist(first.adj)
    second.cor<-second.p<-first.cor*0
    message('Calulate second order partial correlation...')
    if(progressbar)
      pb <- txtProgressBar(min = 0, max =nrow(first2second.edge), style = 3)
    for(i in 1:nrow(first2second.edge)){
      Sys.sleep(0)
      node1<-as.vector(first2second.edge$from[i])
      node2<-as.vector(first2second.edge$to[i])
      sub.edge<-shared.neighbour(node1,node2,first2second.edge,verbose = F)
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
    second.adj<-second.p
    second.adj[second.adj>=p.cut]<-0
    diag(second.adj)<-0
    second.adj[second.adj>0]<-1
    ##
    results<-list(cor=second.cor,pval=second.p,adj=second.adj)
    results.list<-c(results.list,setNames(list(results),'up2second'))
  }
  return(results.list)
}

