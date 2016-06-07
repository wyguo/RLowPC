table.validate<-function(inf.adj,ref.adj,directed=F){
  if(directed){
    table.output<-minet::validate(inf.adj,ref.adj)
    table.zero<-table.output[table.output$thrsh==0,]
    table.output<-rbind(table.output[table.output$thrsh>0,],table.zero[1,])
  } else {
    inf.adj<-pmax(inf.adj,t(inf.adj))
    ref.adj<-pmax(ref.adj,t(ref.adj))
    ref.adj[lower.tri(ref.adj)]<-inf.adj[lower.tri(inf.adj)]<-0
    table.output<-minet::validate(inf.adj,ref.adj)
    tn.remove<-(ncol(ref.adj)^2-ncol(ref.adj))/2+ncol(ref.adj)
    table.output$tn<-table.output$tn-tn.remove
    table.output<-table.output[table.output$tn>=0,]
  }
  return(table.output)
}

