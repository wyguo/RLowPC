#' Relevance low order partial correlation
#'
#'Relevance low order partial correlation (RLowPC) is an improved version of partial correlation (PC)[1, 2].  Instead of removing all remained controls for pair-wise
#'PC calculation, RLowPC selects and regresses the the most relevant controls. See Details.
#'
#' @param data.exp gene expression matrix. Columns are variables and rows are samples.
#' @param edgelist edge list. First column are the name of regulators, second coloumn are the target genes and the third column are the edge weights.
#' @param method a character string to indicate which method is used to calculate correlation. Options are "pearson", "spearman"
#' and "kendall". If \code{pc.estimator="shrink"}, the method is set to "pearson" since there are no "spearman" and "kendall" options for the function \code{\link{pcor.shrink}}
#' to estimate shrink PC.
#' @param pc.estimator a character string to indicate the estimator used to calculate the PC for each pair of nodes in the edge list. Options are "shrink" and "pc",
#' correspoinding to the item (c) and (d) in Details: Step 2, respectively.
#' @param progressbar logical. If TRUE, a progressbar will show to indicate the code runing percentage.
#' @details
#' In general cases, the PC for a pair of genes is calculated by removing all the remained genes (controls). However, there may be a number of inrrelevant controls involved in the removed genes that
#' do not truely connect to the pair of genes. We developed a RLowPC method to calculate PC by only removing more revelant controls. The method is used to refine a pre-inferred network structure by
#' reducing the indirect edges that have been predicted as direct.
#'
#' Step 1: Input a proper size of pre-inferred network, which has room to be improved, for search space of indirect edges. For example the top weighted edges in a inferred PC network can be used as
#' search space. Each pair of genes are assumed to connect to their most relevant neighbour genes since the connecting edges are all top ranked in the whole network.
#'
#' Step 2: Calculate relevance low order partial correlation.  For each pair of genes connected by an edge in the searching
#' space, the edge weight is redefined as (a) zero order PC (correlation) if they do not connect to the same set of neighbour genes, (b) PC by removing all the shared neighbours simultaneously
#' and (c) shrink PC if the covariance matrix used to estimate PC in (b) is not positive definite or invertible. If the searching space is very large, there might still
#' be a number of irrelevant controls involved in shrink PC procedure in (c). An alternative is (d) deleting less connected neighbour genes until the covariance
#' matrix in (b) is positive definite and invertible.
#' @return \code{RLowPC} regurns a new edge list with an additional column of RLowPC edge weights.
#' @references
#' 1. Markowetz F, Spang R: Inferring cellular networks-a review. BMC Bioinformatics 2007, 8 Suppl 6:S5.
#'
#' 2. Sch\"afer J, Strimmer K: A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. Statistical Applications
#' @export
#' @examples
#' ##load library
#' library(RLowPC)
#' library(corpcor)
#' ##load data
#' data(gnwdata)
#' data.exp<-gnwdata$size100$ts1[1:63,]
#' genes<-colnames(data.exp)[-c(1:3)]
#' ##load reference network
#' ref.edge<-gnwdata$size100$net1
#' ref.edge[,3]<-1
#' ref.adj<-edgelist2adjmatrix(edgelist = ref.edge,genes = genes,directed = F)
#'
#'
#' ##filter low expressed genes
#' data2anova<-data.frame(time=factor(paste0(data.exp$experiment,'_',data.exp$time)),data.exp[,-c(1:3)])
#' data.new<-anova2de(data.exp = data2anova,ncol.idx =1,model = 'expression~time',pval.cut = 0.01)
#' data.exp<-data.new$de.ts
#' genes<-data.new$de.gene
#' ref.adj<-ref.adj[genes,genes]
#'
#' ##infer correlation network
#' inf.cor<-abs(cor(data.exp))
#' diag(inf.cor)<-0
#' ##infer PC network
#' inf.pcor<-abs(pcor.shrink(data.exp)[1:length(genes),1:length(genes)])
#' diag(inf.pcor)<-0
#'
#' ##infer LowPC
#' inf.LowPC.edge<-LowPC(data.exp,cutoff = 0.01,cutat='pval')
#' if(is.null(inf.LowPC.edge$secondPC)){
#'   inf.LowPC<-ref.adj*0 } else {
#'     inf.LowPC<-edgelist2adjmatrix(inf.LowPC.edge$secondPC[,1:3],genes = genes,directed = F)
#'   }
#'
#' ##inf RLowPC
#' reduction.sapce<-na.omit(adjmatrix2edgelist(adjmatrix = inf.pcor,directed = F,order = T)[1:1000,])
#' inf.RLowPC.edge<-RLowPC(data.exp = data.exp,edgelist = reduction.sapce,
#'                         method = 'pearson',pc.estimator = 'shrink')
#' inf.RLowPC.edge$cor.weight<-abs(inf.RLowPC.edge$cor.weight)
#' inf.RLowPC<-edgelist2adjmatrix(inf.RLowPC.edge[,c(1,2,4)],genes = genes,directed = T)
#' inf.RLowPC<-abs(inf.RLowPC)
#' inf.RLowPC<-pmax(inf.RLowPC,t(inf.RLowPC))
#'
#' ##infer first order PC based on reduction sapce
#' ###first PC
#' inf.firstPC.edge<-firstPC(data.exp = data.exp,edgelist = reduction.sapce,
#'                           method = 'pearson',controlist = NULL)
#' inf.firstPC<-edgelist2adjmatrix(inf.firstPC.edge[,1:3],genes = genes,directed = F)
#' ###second PC
#' inf.secondPC.edge<-secondPC(data.exp = data.exp,edgelist = reduction.sapce,
#'                             method = 'pearson',controlist = NULL)
#' inf.secondPC<-edgelist2adjmatrix(inf.secondPC.edge[,1:3],genes = genes,directed = F)
#' ##Put the inferred networks into a list.
#' inf.list<-list(Cor=inf.cor,PC=inf.pcor,LowPC=inf.LowPC,firstPC=inf.firstPC,
#'                secondPC=inf.secondPC,RLowPC=inf.RLowPC)
#' sapply(inf.list,dim)
#' dim(ref.adj)
#' ##calculate confusion table
#' inf.table<-lapply(inf.list,function(x) table.evaluate(x,ref.adj = ref.adj))
#' ##cut the tables at top 1000 predictions
#' inf.table<-lapply(inf.table,function(x) table.cut(input.table = x)$table.output)
#'
#' ##plot the figure
#' x11()
#' plotAUC(table.methods = inf.table,
#'         roc.lim = c(0,1,0,1),
#'         pr.lim = c(0,1,0,0.4),
#'         top.precision.lim = c(0,700,0,0.4),
#'         fpr.precision.lim = c(0,1,0,0.4),
#'         color = c('blue','purple','green','lightblue','red','yellow'),
#'         lwd=1,fill = T)




RLowPC<-function(data.exp,edgelist,method='pearson',pc.estimator='shrink',progressbar=T){
  cov2pcor<-function (V)
  {
    ans <- -cov2cor(solve(V))
    diag(ans) <- -diag(ans)
    ans
  }

  if(pc.estimator=='shrink'){
    method='pearson'
  }
  colnames(edgelist)<-c('from','to','weight')
  message(paste0('PC is calculated with method: ',method, ' and pc.estimator: ', pc.estimator,'...'))

  message(paste0('Calculate partial correlation of ', nrow(edgelist), ' pair of nodes ...'))
  t1 <- Sys.time()
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
      cov.matrix<-cov(sub.ts,method = method)
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
          sub.cov<-cov(data.exp[,c(from,to,sub.nodes)],method = method)
          while(!matrixcalc::is.positive.definite(sub.cov) | det(sub.cov) < .Machine$double.eps){
            mean.weight<-mean.weight[-1,]
            sub.nodes<-as.vector(mean.weight[,1])
            sub.cov<-cov(data.exp[,c(from,to,sub.nodes)],method = method)
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
  t2 <- Sys.time()
  message(paste0('Done! Time taken:',round(as.numeric(difftime(t2,t1)),4),' ',units(difftime(t2,t1))))
  edgelist$cor.weight<-cor.weight
  return(edgelist)
}
