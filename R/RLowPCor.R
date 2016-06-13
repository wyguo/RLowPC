#' Relevance low order partial correlation
#'
#'Relevance low order partial correlation (RLowPCor) is an improved version of partial correlation (PC)[1, 2].  Instead of removing all remained controls for pair-wise
#'PC calculation, RLowPCor selects and regresses the the most relevant controls. See Details.
#'
#' @param data.exp gene expression matrix. Columns are variables and rows are samples.
#' @param edgelist edge list. First column are the name of regulators, second coloumn are the target genes and the third column are the edge weights.
#' @param estimator a character string indicating which correlation coefficient (or covariance) is to be computed.Options are "pearson", "spearman"
#' and "kendall". If shrinkage method is used to estimate PC, the estimator is set to "pearson".
#' @param pc.estimator a character string indicating which method is used to estimate the PC of nodes connected to shared neigbours. Options are "shrink" and "pc",
#' correspoinding to the item (c) and (d) in Details: Step 2, respectively.
#' @param progressbar logical. If TRUE, a progressbar will show to indicate the code runing percentage.
#' @details
#' Step 1:  Extract a sparse and scale-free topology from pre-inferred networks as an indirect edge search space for RLowPCor. Unlike a fully connected network,
#' the nodes in the sparse network are assumed to connect to their more relevant neighbours. For example, correlation network can be cut with a range of thresholds
#' until it most fit to scale-free topology.
#'
#' Step 2: Calculate relevance low order partial correlation.  For each pair of nodes connected by an edge in the searching
#' space, the edge weight is redefined as (a) Pearson correlation if they do not connect to the same set of neighbour nodes, (b) PC by removing shared neighbours and
#' (c) shrink PC if the covariance matrix used to estimate PC in (b) is not positive definite or invertible. If the searching space is very large, there might still
#' be a number of irrelevant controls involved in shrink PC procedure (c). An alternative is (d) deleting less connected neighbours of the nodes until the covariance
#' matrix in (b) is positive definite and invertible.
#' @return a network matrix
#' @references
#' 1. Markowetz F, Spang R: Inferring cellular networks-a review. BMC Bioinformatics 2007, 8 Suppl 6:S5.
#'
#' 2. Sch\"afer J, Strimmer K: A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Statistical Applications
#' @export
#' @examples
#' #' Install RLowPCor and load R packages
#' ##Install
#' library(devtools)
#' install_github('wyguo/RLowPCor')
#' ###RLowPCor R package
#' library(RLowPCor)
#' ##
#' library(minet)
#' library(Hmisc)
#' library(ppcor)
#' library(corpcor)
#' library(plyr)
#' library(ggplot2)
#' library(gridExtra)
#'
#'
#' #Load DREAm4 datasets and referecne network
#'
#' data(gnwdata)
#' data.exp<-gnwdata$size100$ts1[,-c(1:2)]
#' genes<-colnames(data.exp)
#' ##reference network edge list
#' ref.edge<-gnwdata$size100$net1
#' ref.edge[,3]<-1
#' ref.adj.raw<-edgelist2adjmatrix(ref.edge,genes)
#' ##refernce network matrix
#' ref.adj<-ref.adj.raw[genes,genes]
#'
#'
#'
#'
#' #Infer correlation network
#'
#' inf.cor<-abs(cor(data.exp))
#' diag(inf.cor)<-0
#'
#'
#' #Deterimine to top weighted search space.
#' #The reduction of search space of correlation network is selected by testing the scale-free fitness of top weighted edges the number of which ranging from 100 to 2000. The testing results shows top edge number=180 indicates the optimal case.
#'
#' fitness<-fit.structure(net=inf.cor,top.vector=seq(100,2000,by=20),method='scale-free')
#' top=as.numeric(names(fitness[fitness==max(fitness)]))
#' ##choose the top 180 edges.
#' inf.cor.top<-top.network(inf.cor,top=top)
#'
#'
#' #Infer shrink PC network
#'
#' inf.pcor<-abs(pcor.shrink(data.exp,verbose = F))[1:ncol(data.exp),1:ncol(data.exp)]
#' diag(inf.pcor)<-0
#' ##choose the top 180 edges
#' inf.pcor.top<-top.network(inf.pcor,top=top)
#'
#'
#' #Infer firstPCor and secondPCor
#'
#' ##' use the top weighted structure of correlation network as search space
#' adjmatrix<-inf.cor.top
#' ##' infer first order partial correlation network
#' inf.firstpcor<-abs(firstPCor(adjmatrix = adjmatrix,data.exp,progressbar = F)$cor)
#' ##'  infer second order partial correlation network
#' inf.secondpcor<-abs(secondPCor(adjmatrix = adjmatrix,data.exp,progressbar = F)$cor)
#'
#'
#' #Infer RLowPCor
#'
#' ##RLowPCor
#' inf.edge<-RLowPCor(data.exp,edgelist=adjmatrix2edgelist(inf.cor.top,directed=T),progressbar = F)
#' inf.rlowpcor<-edgelist2adjmatrix(inf.edge[,c(1,2,4)],genes=colnames(inf.cor.top),directed=F)
#'
#'
#' #Infer LowPCor
#'
#' inf.lopcor<-LowPCor(data.exp,p.cut = 0.001,progressbar = F)$up2second$cor
#' inf.lopcor.top<-top.network(inf.lopcor,top=top)
#'
#' #Generate confusion table
#'
#' ##convert the inferred metwork matrix to symmetric to match with reference network
#' inf.cor.top<-pmax(inf.cor.top,t(inf.cor.top))
#' inf.pcor.top<-pmax(inf.pcor.top,t(inf.pcor.top))
#' inf.firstpcor<-pmax(inf.firstpcor,t(inf.firstpcor))
#' inf.secondpcor<-pmax(inf.secondpcor,t(inf.secondpcor))
#' inf.lopcor.top<-pmax(inf.lopcor.top,t(inf.lopcor.top))
#' ##generate confusion table
#' table.firstpcor<-table.evaluate(inf.firstpcor,ref.adj)
#' table.secondpcor<-table.evaluate(inf.secondpcor,ref.adj)
#' table.cor<-table.evaluate(inf.cor.top,ref.adj)
#' table.pcor<-table.evaluate(inf.pcor.top,ref.adj)
#' table.rlowpcor<-table.evaluate(inf.rlowpcor,ref.adj)
#' table.lowpcor<-table.evaluate(inf.lopcor.top,ref.adj)
#'
#'
#' #plot the partial auroc and aupr curves
#'
#' cat(paste0('\n Predictions of top ',top,' edges\n'))
#' f=plotAUC(list(Cor=table.cor,
#'                firstPCor=table.firstpcor,
#'                secondPCor=table.secondpcor,
#'                shrink.PCor=table.pcor,
#'                RLowPCor=table.rlowpcor,
#'                LowPCor=table.lowpcor),lwd=1,fill=T,pr.lim = c(0,0.25,0,0.5),roc.lim = c(0,0.04,0,0.25),color = rainbow(6))



RLowPCor<-function(data.exp,edgelist,estimator='pearson',pc.estimator='shrink',progressbar=T){
  cov2pcor<-function (V)
  {
    ans <- -cov2cor(solve(V))
    diag(ans) <- -diag(ans)
    ans
  }

  if(pc.estimator=='shrink'){
    estimator='pearson'
  }
  colnames(edgelist)<-c('from','to','weight')
  message(paste0('PC is calculated with estimator: ',estimator, ' and pc.estimator: ', pc.estimator,'...'))

  message(paste0('Calculate partial correlation of ', nrow(edgelist), ' pair of nodes ...'))
  if(progressbar)
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
      sub.ts<-data.exp[,sub.genes]
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
    if(progressbar)
      setTxtProgressBar(pb,i)
  }
  if(progressbar)
    close(pb)
  edgelist$cor.weight<-cor.weight
  return(edgelist)
}
