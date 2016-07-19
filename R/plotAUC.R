#' Plot AUC curves
#'
#' The function is used to plot AUROC, AUPR, pAUROC and pAUPR curves.
#' @param table.methods a list of confusion tables of comparison for multiple networks to the same reference network. The tables are the outputs of \code{\link{table.evalution}}.
#' The names of the list \code{names(table.methods)} will be used
#' to plot the legend of the curves.
#' @param plot.method a character string indicates to plot AUROC ("auroc"), AUPR ("aupr"), precision vs top weighted edges ("top.precision") or
#' precision vs false positive rate ("fpr.precision") curves.
#' @param roc.lim a vector to specify the plot limits, i.e. xmin, xmax, ymin and ymax of AUROC curves.
#' @param pr.lim a vector to specify the plot limits, i.e. xmin, xmax, ymin and ymax of AUPR curves.
#' @param top.precision.lim a vector to specify the plot limits, i.e. xmin, xmax, ymin and ymax of precision vs top weighted edges curves.
#' @param fpr.precision.lim a vector to specify the plot limits, i.e. xmin, xmax, ymin and ymax of precision vs false positive rate curves.
#' @param fill logical. If TRUE, the area under the curvers will be filled with colors.
#' @param color a color string vectors to specify the colors for the curves.
#' @param plot.ncol the number of columns of the layout of mutliple plots.
#' @param ... additional options for \code{geom_line(...)}.
#' @examples
#' ##load librarys
#' library(RLowPC)
#' library(minet)
#' library(gridExtra)
#' library(ggplot2)
#' ##get dream4 datasets
#' data(gnwdata)
#' data.exp<-gnwdata$size100$ts1[,-c(1:3)]
#' genes<-colnames(data.exp)
#' ref.edge<-gnwdata$size100$net1
#' ref.edge[,3]<-1
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
#' plotAUC(table.methods,fill=T,lwd=1)
#'
#' @return \code{plotAUC} returns plots.
#' @export

plotAUC<-function(table.methods,plot.method=c('auroc','aupr','top.precision','fpr.precision'),
                  roc.lim=c(0,1,0,1),
                  pr.lim=c(0,1,0,1),
                  top.precision.lim=c(0,1000,0,1),
                  fpr.precision.lim=c(0,1,0,1),
                  fill=F,
                  color=NA,
                  plot.ncol=2,...){
  library(ggplot2)
  get.legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  top.precision.frame<-fpr.precision.frame<-roc.frame<-pr.frame<-data.frame()
  for(i in 1:length(table.methods)){
    table.conf<-confusion(table.methods[[i]])
    table.conf<-na.omit(table.conf)
    fpr<-table.conf$fpr
    tpr<-table.conf$tpr
    recall<-table.conf$recall
    precision<-table.conf$precision
    top<-table.conf$tp+table.conf$fp
    top.precision<-table.conf$precision
    fpr.fpr<-table.conf$fpr
    fpr.precision<-table.conf$precision
    if(fill){
      fpr<-c(0,fpr,max(fpr),max(fpr))
      tpr<-c(0,tpr,0,0)
      recall<-c(0,recall,max(recall),0)
      precision<-c(0,precision,0,0)
      top<-c(top.precision.lim[1],top,max(top),0)
      top.precision<-c(0,top.precision,0,0)
      fpr.fpr<-c(fpr.precision.lim[1],fpr.fpr,max(fpr.fpr),0)
      fpr.precision<-c(0,fpr.precision,0,0)
    }
    roc.frame<-rbind(roc.frame,data.frame(network=names(table.methods)[i],fpr=fpr,tpr=tpr))
    pr.frame<-rbind(pr.frame,data.frame(network=names(table.methods)[i],recall=recall,precision=precision))
    top.precision.frame<-rbind(top.precision.frame,
                               data.frame(network=names(table.methods)[i],
                                          top=top,
                                          precision=top.precision
                               )
    )
    fpr.precision.frame<-rbind(fpr.precision.frame,
                               data.frame(network=names(table.methods)[i],
                                          fpr=fpr.fpr,
                                          precision=fpr.precision)
    )
  }
  plot.f<-list()

  ###Title of the plots
  title.roc<-'AUROC'
  title.pr<-'AUPR'
  if(roc.lim[2]<1)
    title.roc<-paste0('p',title.roc)
  if(pr.lim[2]<1)
    title.pr<-paste0('p',title.pr)
  ###auroc
  if('auroc' %in% plot.method){
    random.data<-data.frame(network='random',fpr=c(roc.lim[1],roc.lim[2]),tpr=c(roc.lim[1],roc.lim[2]))
    f.roc<-ggplot(aes(x=fpr,y=tpr,color=network,group=network),data=roc.frame)+geom_line(...)+theme_bw()+
      labs(x='FPR',y='TPR',title=title.roc)+
      coord_cartesian(xlim=c(0,roc.lim[2]),ylim=c(0,roc.lim[4]))+
      theme(legend.position='bottom')+
      guides(fill=guide_legend(nrow=2,byrow=F))+
      geom_line(aes(x=fpr,y=tpr),data = random.data,color='black',linetype='longdash',...)+
      geom_text(x=0.9*roc.lim[2],y=0.75*roc.lim[2],label='random',color='black')
    ##color
    if(all(!is.na(color))){
      f.roc<-f.roc+scale_colour_manual(values = color)+
        scale_fill_manual(values = color)
    }


    leg<-get.legend(f.roc)
    f.roc<-f.roc+theme(legend.position='none')
    if(fill){
      f.roc<-f.roc+geom_polygon(data=roc.frame, mapping=aes(x=fpr, y=tpr, group=network,fill=network),alpha=0.05)
    }
    plot.f<-c(plot.f,auroc=list(f.roc))
  }


  ###aupr
  if('aupr' %in% plot.method){
    f.pr<-ggplot(aes(x=recall,y=precision,color=network,group=network),data=pr.frame)+geom_line(...)+theme_bw()+
      labs(x='Recall',y='Precision',title=title.pr)+
      coord_cartesian(xlim=c(0,pr.lim[2]),ylim=c(0,pr.lim[4]))+
      theme(legend.position='bottom')+
      guides(fill=guide_legend(nrow=2,byrow=F))
    ##color
    if(all(!is.na(color))){
      f.pr<-f.pr+scale_colour_manual(values = color)+
        scale_fill_manual(values = color)
    }
    leg<-get.legend(f.pr)
    f.pr<-f.pr+theme(legend.position='none')
    if(fill){
      f.pr<-f.pr+geom_polygon(data=pr.frame, mapping=aes(x=recall, y=precision, group=network,fill=network),alpha=0.05)
    }
    plot.f<-c(plot.f,aupr=list(f.pr))
  }

  ###top vs precision
  if('top.precision' %in% plot.method){
    f.top.precision<-ggplot(aes(x=top,y=precision,color=network,group=network),data=top.precision.frame)+geom_line(...)+theme_bw()+
      labs(x='Top weighted edges',y='Precision',title='Precision at top edge cutoffs')+
      coord_cartesian(xlim=c(0,top.precision.lim[2]),ylim=c(0,top.precision.lim[4]))+
      theme(legend.position='bottom')+
      guides(fill=guide_legend(nrow=2,byrow=F))
    ##color
    if(all(!is.na(color))){
      f.top.precision<-f.top.precision+scale_colour_manual(values = color)+
        scale_fill_manual(values = color)
    }
    leg<-get.legend(f.top.precision)
    f.top.precision<-f.top.precision+theme(legend.position='none')
    if(fill){
      f.top.precision<-f.top.precision+geom_polygon(data=top.precision.frame, mapping=aes(x=top, y=precision, group=network,fill=network),alpha=0.05)
    }
    plot.f<-c(plot.f,aupr=list(f.top.precision))
  }


  ###fpr vs precision
  if('fpr.precision' %in% plot.method){
    f.fpr.precision<-ggplot(aes(x=fpr,y=precision,color=network,group=network),data=fpr.precision.frame)+geom_line(...)+theme_bw()+
      labs(x='FPR',y='Precision',title='Precision at FPR cutoffs')+
      coord_cartesian(xlim=c(0,fpr.precision.lim[2]),ylim=c(0,fpr.precision.lim[4]))+
      theme(legend.position='bottom')+
      guides(fill=guide_legend(nrow=2,byrow=F))
    ##color
    if(all(!is.na(color))){
      f.fpr.precision<-f.fpr.precision+scale_colour_manual(values = color)+
        scale_fill_manual(values = color)
    }
    leg<-get.legend(f.fpr.precision)
    f.fpr.precision<-f.fpr.precision+theme(legend.position='none')
    if(fill){
      f.fpr.precision<-f.fpr.precision+geom_polygon(data=fpr.precision.frame, mapping=aes(x=fpr, y=precision, group=network,fill=network),alpha=0.05)
    }
    plot.f<-c(plot.f,aupr=list(f.fpr.precision))
  }


  plots<-gridExtra::grid.arrange(do.call(gridExtra::arrangeGrob,c(plot.f,ncol=plot.ncol)),nrow=2,leg,heights=c(10,2))
  #plots<-gridExtra::grid.arrange(gridExtra::marrangeGrob(plot.f,ncol=2,nrow=3),leg)
  message('Plots done!')
  return(plots)
}
