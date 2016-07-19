# RLowPC
R package to construct Relevance Low order Partial Correlation gene networks
##Installation
To install the package, run the R codes
```{r}
library(devtools)
install_github('wyguo/RLowPC')
```
Additional packages used

ppcor,corpcor, minet, netbenchmark, ggplot2, gridExtra, plyr, fdrtool

##Examples of RLowPC R package
####Example1: RLowPC function
```{r,cache=T,include=T,message='hide',warning=F,eval=F}
##load library
library(RLowPC)
library(corpcor)
##load data
data(gnwdata)
data.exp<-gnwdata$size100$ts1[1:63,]
genes<-colnames(data.exp)[-c(1:3)]
##load reference network
ref.edge<-gnwdata$size100$net1
ref.edge[,3]<-1
ref.adj<-edgelist2adjmatrix(edgelist = ref.edge,genes = genes,directed = F)


##filter low expressed genes
data2anova<-data.frame(time=factor(paste0(data.exp$experiment,'_',data.exp$time)),data.exp[,-c(1:3)])
data.new<-anova2de(data.exp = data2anova,ncol.idx =1,model = 'expression~time',pval.cut = 0.01)
data.exp<-data.new$de.ts
genes<-data.new$de.gene
ref.adj<-ref.adj[genes,genes]

##infer correlation network
inf.cor<-abs(cor(data.exp))
diag(inf.cor)<-0
##infer PC network
inf.pcor<-abs(pcor.shrink(data.exp)[1:length(genes),1:length(genes)])
diag(inf.pcor)<-0

##infer LowPC
inf.LowPC.edge<-LowPC(data.exp,cutoff = 0.01,cutat='pval')
if(is.null(inf.LowPC.edge$secondPC)){
  inf.LowPC<-ref.adj*0 } else {
    inf.LowPC<-edgelist2adjmatrix(inf.LowPC.edge$secondPC[,1:3],genes = genes,directed = F)
  }

##inf RLowPC
reduction.sapce<-na.omit(adjmatrix2edgelist(adjmatrix = inf.pcor,directed = F,order = T)[1:1000,])
inf.RLowPC.edge<-RLowPC(data.exp = data.exp,edgelist = reduction.sapce,
                        method = 'pearson',pc.estimator = 'shrink')
inf.RLowPC.edge$cor.weight<-abs(inf.RLowPC.edge$cor.weight)
inf.RLowPC<-edgelist2adjmatrix(inf.RLowPC.edge[,c(1,2,4)],genes = genes,directed = T)
inf.RLowPC<-abs(inf.RLowPC)
inf.RLowPC<-pmax(inf.RLowPC,t(inf.RLowPC))

##infer first order PC based on reduction sapce
###first PC
inf.firstPC.edge<-firstPC(data.exp = data.exp,edgelist = reduction.sapce,
                          method = 'pearson',controlist = NULL)
inf.firstPC<-edgelist2adjmatrix(inf.firstPC.edge[,1:3],genes = genes,directed = F)
###second PC
inf.secondPC.edge<-secondPC(data.exp = data.exp,edgelist = reduction.sapce,
                            method = 'pearson',controlist = NULL)
inf.secondPC<-edgelist2adjmatrix(inf.secondPC.edge[,1:3],genes = genes,directed = F)
##Put the inferred networks into a list.
inf.list<-list(Cor=inf.cor,PC=inf.pcor,LowPC=inf.LowPC,firstPC=inf.firstPC,
               secondPC=inf.secondPC,RLowPC=inf.RLowPC)
sapply(inf.list,dim)
dim(ref.adj)
##calculate confusion table
inf.table<-lapply(inf.list,function(x) table.evaluate(x,ref.adj = ref.adj))
##cut the tables at top 1000 predictions
inf.table<-lapply(inf.table,function(x) table.cut(input.table = x)$table.output)

##plot the figure
x11()
plotAUC(table.methods = inf.table,
        roc.lim = c(0,1,0,1),
        pr.lim = c(0,1,0,0.4),
        top.precision.lim = c(0,700,0,0.4),
        fpr.precision.lim = c(0,1,0,0.4),
        color = c('blue','purple','green','lightblue','red','yellow'),
        lwd=1,fill = T)


```

####Example2: average.consensus function
```{r,cache=T,include=T,message='hide',warning=F}
 ##create two random networks
 library(RLowPC)
 set.seed(4)
 net1<-abs(matrix(rnorm(16),4,4))
 net1<-pmax(net1,t(net1))
 diag(net1)<-0
 set.seed(5)
 net2<-abs(matrix(rnorm(16),4,4))
 net2<-pmax(net2,t(net2))
 diag(net2)<-0
 dimnames(net1)<-dimnames(net2)<-list(letters[1:4],letters[1:4])
 ##put networks into list
 net.list<-list(net1=net1,net2=net2)
 inf.consensus<-average.consensus(adjmatrix.list = net.list,directed = F)
 adj2rankadj(net1)
 adj2rankadj(net2)
 inf.consensus
```

####Example3: table.evaluate function
```{r,cache=T,include=T,message='hide',warning=F}
##load library
library(RLowPC)
##load data
data(gnwdata)
data.exp<-gnwdata$size100$ts1[,-c(1:3)]
genes<-colnames(data.exp)
ref.edge<-gnwdata$size100$net1
ref.edge[,3]<-1
ref.adj<-edgelist2adjmatrix(ref.edge,genes)
inf.cor<-abs(cor(data.exp))
diag(inf.cor)<-0
table.cor<-table.evaluate(inf.adj = inf.cor,ref.adj = ref.adj)
head(table.cor)
 
```
####Example4: edgelist2adjmatrix function
```{r,cache=T,include=T,message='hide',warning=F}
##load data
library(RLowPC)
data(gnwdata)
data.exp<-gnwdata$size100$ts1[,-c(1:3)]
genes<-colnames(data.exp)
ref.edge<-gnwdata$size100$net1
ref.edge[,3]<-1
ref.adj<-edgelist2adjmatrix(ref.edge,genes,directed=F)
 
```

####Example5: adjmatrix2edgelist function
```{r,cache=T,include=T,message='hide',warning=F}
##load data
library(RLowPC)
data(gnwdata)
data.exp<-gnwdata$size100$ts1[,-c(1:3)]
genes<-colnames(data.exp)
##build correlation network
inf.cor<-abs(cor(data.exp))
diag(inf.cor)<-0
##convert matrix to edge list
adjmatrix2edgelist(inf.cor)
 
```

