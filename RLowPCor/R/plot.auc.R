#' Plot AUC curves
#'
#' The function is used to plot AUROC, AUPR, pAUROC and pAUPR curves.
#' @param table.methods a confusion table list of multiple methods, which are from the output of \code{\link{table.evalution}}. The names of the list will be used
#' to plot the legend of the curves.
#' @param plot.method a character string indicates to plot "auroc" or "aupr" curves.
#' @param roc.lim a vector to tell the plot limits, i.e. xmin,xmax,ymin,ymax, of AUROC curves.
#' @param pr.lim a vector to tell the plot limits, i.e. xmin,xmax,ymin,ymax, of AUPR curves.
#' @param fill logical. If TRUE, the area under the curvers will be filled with colors.
#' @details The fpr or recall cutoffs to make pAUROC or pAUPR curves can be generate from the funciton \code{\link{confusion}} and specify in \code{roc.lim} and \code{pr.lim}
#' @examples
#' ##load librarys
#' library(networkBMA)
#' library(RLowPCor)
#' library(minet)
#' library(gridEXtra)
#' library(gglot2)
#' ##get dream4 datasets
#' data(dream4)
#' data.exp<-dream4ts100[[1]][,-c(1:2)]
#' genes<-colnames(data.exp)
#' ref.edge<-dream4gold100[[1]]
#' ref.adj<-edgelist2adjmatrix(ref.edge,genes)
#' ##infer gene networks
#' inf.cor<-abs(cor(data.exp))
#' diag(inf.cor)<-0
#' inf.mi<-build.mim(data.exp)
#' inf.clr<-clr(inf.mi)
#' inf.mrnet<-mrnet(inf.mi)
#' ##generate confusion tables
#' table.cor<-table.evaluate(inf.adj = inf.cor,ref.adj = ref.adj)
#' table.mi<-table.evaluate(inf.adj=inf.mi,ref.adj = ref.adj)
#' table.clr<-table.evaluate(inf.adj=inf.clr,ref.adj = ref.adj)
#' table.mrnet<-table.evaluate(inf.adj=inf.mrnet,ref.adj = ref.adj)
#' ##put confusion tables into list, and set names as the methods
#' table.methods<-list(cor=table.cor,mi=table.mi,clr=table.clr,mrnet=table.mrnet)
#' plot.auc(table.methods,plot.method = 'auroc',roc.lim = c(0,0.5,0,1))
#'
#' @return Plots.
#' @export

plot.auc<-function(table.methods,plot.method=c('auroc','aupr'),roc.lim=c(0,1,0,1),pr.lim=c(0,1,0,1),fill=F){
  get.legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  roc.frame<-pr.frame<-data.frame()
  for(i in 1:length(table.methods)){
    table.conf<-confusion(table.methods[[i]])
    fpr<-table.conf$fpr
    tpr<-table.conf$tpr
    recall<-table.conf$recall
    precision<-table.conf$precision
    if(fill){
      fpr<-c(roc.lim[2],fpr,roc.lim[2])
      tpr<-c(0,tpr,0)
      recall<-c(pr.lim[2],recall,pr.lim[2])
      precision<-c(0,precision,0)

    }
    roc.frame<-rbind(roc.frame,data.frame(network=names(table.methods)[i],fpr=fpr,tpr=tpr))
    pr.frame<-rbind(pr.frame,data.frame(network=names(table.methods)[i],recall=recall,precision=precision))
  }
  plot.f<-list()
  if('auroc' %in% plot.method){
    f.roc<-ggplot(aes(x=fpr,y=tpr,color=network,group=network),data=roc.frame)+geom_line(lwd=1)+theme_bw()+
      geom_segment(x = 0, y = 0,xend=roc.lim[2],yend=roc.lim[2],lwd=1,linetype = "longdash",col='gray5',stat = 'identity')+
      labs(x='FPR',y='TPR',title='AUROC')+
      xlim(c(roc.lim[1],roc.lim[2]))+
      ylim(c(roc.lim[3],roc.lim[4]))+
      theme(legend.position='bottom')+
      guides(fill=guide_legend(nrow=2,byrow=F))
    leg<-get.legend(f.roc)
    f.roc<-f.roc+theme(legend.position='none')
    if(fill){
      f.roc<-f.roc+geom_polygon(data=roc.frame, mapping=aes(x=fpr, y=tpr, group=network,fill=network),alpha=0.1)
    }
    plot.f<-c(plot.f,auroc=list(f.roc))
  }
  if('aupr' %in% plot.method){
    f.pr<-ggplot(aes(x=recall,y=precision,color=network,group=network),data=pr.frame)+geom_line(lwd=1)+theme_bw()+
      labs(x='Recall',y='Precision',title='AUPR')+
      xlim(c(pr.lim[1],pr.lim[2]))+
      ylim(c(pr.lim[3],pr.lim[4]))+
      theme(legend.position='bottom')+
      guides(fill=guide_legend(nrow=2,byrow=F))
    leg<-get.legend(f.pr)
    f.pr<-f.pr+theme(legend.position='none')
    if(fill){
      f.pr<-f.pr+geom_polygon(data=pr.frame, mapping=aes(x=recall, y=precision, group=network,fill=network),alpha=0.1)
    }
    plot.f<-c(plot.f,aupr=list(f.pr))
  }
  gridExtra::grid.arrange(do.call(gridExtra::arrangeGrob,c(plot.f,nrow=1)),leg,heights=c(10, 1))
  message('Plots done!')
  return()
}
