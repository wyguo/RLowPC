#' Relevance low order partial correlation
#'
#' Consensus network is built of taking the average ranks of the edges from multiple network predictions.
#'
#' @param data.exp gene expression matrix. Columns are variables and rows are samples.
#' @param edgelist edge list. First column are the name of regulators, second coloumn are the target genes and the third column are the edge weights.
#' @param estimator a character string indicating which correlation coefficient (or covariance) is to be computed.Options are "pearson", "spearman"
#' and "kendall". If shrinkage method is used to estimate PC, the estimator is set to "pearson".
#' @param pc.estimator a character string indicating which method is used to estimate the PC of nodes connected to shared neigbours. Options are "shrink" and "pc",
#' correspoinding to the item (c) and (d) in Details: Step 2, respectively.
#' @details
#' Step 1:  Extract a sparse and scale-free topology from pre-inferred networks as an indirect edge search space for RLowPCor. Unlike a fully connected network,
#' the nodes in the sparse network are assumed to connect to their more relevant neighbours. For example, correlation network can be cut with a range of thresholds
#' until it most fit to scale-free topology. Step 2: Calculate relevance low order partial correlation.  For each pair of nodes connected by an edge in the searching
#' space, the edge weight is redefined as (a) Pearson correlation if they do not connect to the same set of neighbour nodes, (b) PC by removing shared neighbours and
#' (c) shrink PC if the covariance matrix used to estimate PC in (b) is not positive definite or invertible. If the searching space is very large, there might still
#' be a number of irrelevant controls involved in shrink PC procedure (c). An alternative is (d) deleting less connected neighbours of the nodes until the covariance
#' matrix in (b) is positive definite and invertible.
#' @return A network with rank weighted edges.
#' @references Sch\"afer J, Strimmer K: A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Statistical Applications
#' in Genetics and Molecular Biology, The Berkeley Electronic Press 2005, 4(1).
#' @export
#' @examples
#' ##load size 100_1 network DREAM4 datasets
#' library(networkBMA)
#' library(RLowPCor)
#' data(dream4)
#' data.exp<-dream4ts100[[1]]
#' #create edge list
#' edgelist<-dream4gold100[[1]]
#' edgelist<-edgelist[edgelist[,3]>0,]
#' ##infer RLowPCor network
#' inf.net=RLowPCor(data.exp = data.exp[,-c(1:2)],edgelist = edgelist)

RLowPCor<-function(data.exp,edgelist,estimator='pearson',pc.estimator='shrink'){
  cov2pcor<-function (V)
  {
    ans <- -cov2cor(solve(V))
    diag(ans) <- -diag(ans)
    ans
  }

  if(pc.estimator=='shrink'){
    estimator='pearson'
    cat('\n The shrink PC is estimated with Pearson algorithm. \n')
  }
  colnames(edgelist)<-c('from','to','weight')
  if('shrink' %in% pc.estimator){
    cat('\n PC is calucated by shrink method \n')
  }
  if('pc' %in% pc.estimator){
    cat('\n PC is calcuated using  general method \n')
  }

  cat('\n Calculate partial correlation of ', nrow(edgelist), 'pair of nodes ...')
  pb <- txtProgressBar(min = 0, max =nrow(edgelist), style = 3)
  cor.weight<-vector()
  for(i in 1:nrow(edgelist)){
    Sys.sleep(0)
    from=as.vector(edgelist$from[i])
    to=as.vector(edgelist$to[i])
    sub.edgelist<-edgelist[which(edgelist$to %in% from | edgelist$to %in% to),]
    filter.gene<-table(sub.edgelist[,1])
    filter.gene<-names(filter.gene)[filter.gene==1]
    if(length(filter.gene)==nrow(sub.edgelist)){
      cor.weight<-c(cor.weight,cor(data.exp[,from],data.exp[,to]))
    } else {
      if(length(filter.gene)==0) {
        sub.edgelist.new<-sub.edgelist
      } else {
        sub.edgelist.new<-sub.edgelist[-which(sub.edgelist[,1] %in% filter.gene),]
      }
      sub.genes<-unique(c(as.vector(sub.edgelist.new[,1]),as.vector(sub.edgelist.new[,2])))
      sub.ts<-data.exp[,sub.genes][1:210,]
      cov.matrix<-cov(sub.ts,method = estimator)
      if(matrixcalc::is.positive.definite(cov.matrix) & det(cov.matrix) >= .Machine$double.eps){
        inf.pcor<-cov2pcor(cov.matrix)
      } else {
        if(pc.estimator=='shrink'){
          inf.pcor<-pcor.shrink(sub.ts,verbose = F)[1:ncol(sub.ts),1:ncol(sub.ts)]
        }
        if(pc.estimator=='pc'){
          mean.weight<-aggregate(sub.edgelist.new$weight,by=list(sub.edgelist.new$from),FUN=mean)
          mean.weight<-mean.weight[order(mean.weight[,2],decreasing = F),]
          sub.nodes<-as.vector(mean.weight[,1])
          sub.cov<-cov(data.exp[,c(from,to,sub.nodes)],method = estimator)
          while(!matrixcalc::is.positive.definite(sub.cov) | det(sub.cov) < .Machine$double.eps){
            mean.weight<-mean.weight[-1,]
            sub.nodes<-as.vector(mean.weight[,1])
            sub.cov<-cov(data.exp[,c(from,to,sub.nodes)],method = estimator)
          }
          inf.pcor<-cov2pcor(sub.cov)
        }
      }
      cor.weight<-c(cor.weight,inf.pcor[from,to])
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)
  edgelist$cor.weight<-cor.weight
  return(edgelist)
}
