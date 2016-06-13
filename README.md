# RLowPCor
R package to construct Relevance Low order Partial Correlation gene networks
##Installation
To install the package, run the R codes
```{r}
library(devtools)
install_github('wyguo/RLowPCor')
```
Additional packages used

ppcor,corpcor, minet, netbenchmark, ggplot2, gridExtra, plyr, Hmisc

##Examples of RLowPCor R package
####Example1: RLowPCor function
```{r,cache=T,include=T,message='hide',warning=F,eval=F}
##load size 100_1 network DREAM4 datasets
 library(networkBMA)
 library(RLowPCor)
 data(dream4)
 data.exp<-dream4ts100[[1]]
 #create edge list
 edgelist<-dream4gold100[[1]]
 edgelist<-edgelist[edgelist[,3]>0,]
 ##infer RLowPCor network
 inf.net=RLowPCor(data.exp = data.exp[,-c(1:2)],edgelist = edgelist)

```

####Example2: average.consensus function
```{r,cache=T,include=T,message='hide',warning=F}
 ##create two random networks
 library(RLowPCor)
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
 ##create two random networks
 library(networkBMA)
 library(RLowPCor)
 data(dream4)
 data.exp<-dream4ts10[[1]][,-c(1:2)]
 genes<-colnames(data.exp)
 ref.edge<-dream4gold10[[1]]
 ref.adj<-edgelist2adjmatrix(ref.edge,genes)
 inf.cor<-abs(cor(data.exp))
 diag(inf.cor)<-0
 table.cor<-table.evaluate(inf.adj = inf.cor,ref.adj = ref.adj)
 head(table.cor)
 
```
####Example4: edgelist2adjmatrix function
```{r,cache=T,include=T,message='hide',warning=F}
 ##create two random networks
 library(networkBMA)
 library(RLowPCor)
 data(dream4)
 data.exp<-dream4ts10[[1]][,-c(1:2)]
 genes<-colnames(data.exp)
 ref.edge<-dream4gold10[[1]]
 ref.adj<-edgelist2adjmatrix(ref.edge,genes)
 ref.adj
 
```


####Example5: adjmatrix2edgelist function
```{r,cache=T,include=T,message='hide',warning=F}
 ##create two random networks
 library(networkBMA)
 library(RLowPCor)
 data(dream4)
 data.exp<-dream4ts10[[1]][,-c(1:2)]
 genes<-colnames(data.exp)
 inf.cor<-abs(cor(data.exp))
 diag(inf.cor)<-0
 adjmatrix2edgelist(inf.cor)
 
```

####Example6: plot.auc function
```{r,cache=T,include=T,message='hide',warning=F}
 library(networkBMA)
 library(RLowPCor)
 library(minet)
 data(dream4)
 data.exp<-dream4ts100[[1]][,-c(1:2)]
 genes<-colnames(data.exp)
 ref.edge<-dream4gold100[[1]]
 ref.adj<-edgelist2adjmatrix(ref.edge,genes)
 inf.cor<-abs(cor(data.exp))
 diag(inf.cor)<-0
 inf.mi<-build.mim(data.exp)
 inf.clr<-clr(inf.mi)
 inf.mrnet<-mrnet(inf.mi)
 table.cor<-table.evaluate(inf.adj = inf.cor,ref.adj = ref.adj)
 table.mi<-table.evaluate(inf.adj=inf.mi,ref.adj = ref.adj)
 table.clr<-table.evaluate(inf.adj=inf.clr,ref.adj = ref.adj)
 table.mrnet<-table.evaluate(inf.adj=inf.mrnet,ref.adj = ref.adj)
 table.methods<-list(cor=table.cor,mi=table.mi,clr=table.clr,mrnet=table.mrnet)
```
