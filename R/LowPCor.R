
#' Calculate low order partial correlation
#'
#' Calculate up to second order partial correlation (PC) as shown in paper [1]. A pair of nodes are connected if the
#' zero order (general correlation), all first order and all second order PC are significantly different to zero.
#'
#' @param data.exp gene expression data matrix with genes in columns and samples in rows.
#' @param p.cut p-value significance threshold for zero, first and second PC.
#' @param up2 a numeric value to indicate the maximum order for recusive PC calculation. If \code{up2=i, (i=0,1,2)} the ouput results are i order
#' PC matrix and corresponding p-value matrix.
#' @param p.quant.prob a probability value for pair-wise p-value quantile to choose as the significance metric of different order PC. For example,
#' in a \eqn{n} size network there are \eqn{n-2} p-values of first order PC for each pair of nodes. The \code{quantile(q.values,probs=p.quant.prob)} of p-values
#' is used as the threshold of signififcane. \code{p.quant.prob=1} present the maximum value of \eqn{n-2} p-values is selected.
#' @param method a string character of method used to estimate correlation. Options are "pearson" and "spearman".
#' @param progressbar logical. If TRUE, a progressbar will show to indicate the code runing percentage.
#' @references
#' 1. Zuo Y, Yu G, Tadesse MG, Ressom HW: Biological network inference using low order partial correlation. Methods (San Diego, Calif) 2014, 69(3):266-273.
#' @return A list contains zero, first or secnd order PC matrix, the corresponding p-value matrix and the adjacency matrix after applying \code{p.cut} to the
#' p-value matrix.
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

