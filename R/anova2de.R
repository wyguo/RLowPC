#' Use ANOVA to filter low epxressed genes
#'
#' Analysis of variance is used to filter the low expressed genes. More details can be seen in \code{\link{aov}}.
#' @param data.exp gene expression data matrix with variables in columns and samples in rows. The first columns are the indeces for comparison, such as
#' time points and conditions for comparision.
#' @param ncol.idx the number of index columns.
#' @param pval.cut a numeric value for significance cut-off. The default is \code{pval.cut=0.05}.
#' @param model the model used for comparison. The default is \code{model="expression~time"}.
#'
#' @return \code{anova2de} returns a list includes the ANOVA p-value results, the names of differentially expressed genes and the time-series data of differential expressed genes.
#' @seealso \code{\link{aov}}, \code{\link{RLowPC}}
#' @export


anova2de<-function(data.exp,ncol.idx,pval.cut=0.05,model='expression~time'){
  gene.names<-colnames(data.exp)[-c(1:ncol.idx)]
  results=data.frame()
  message('ANOVA analysis of differentially expressed genes. Gene number=',length(gene.names))
  pb <- txtProgressBar(min = 0, max =length(gene.names), style = 3)
  for(i in 1:length(gene.names)){
    Sys.sleep(0)
    onegene=data.frame(data.exp[,1:ncol.idx],data.exp[,colnames(data.exp)==gene.names[i]])
    colnames(onegene)<-c(colnames(data.exp)[1:ncol.idx],strsplit(model,'~')[[1]][1])
    fit<-aov(as.formula(model),data=onegene)
    pval=as.numeric(format(summary(fit)[[1]]$"Pr(>F)"[1],digits=4))
    sub.result=data.frame(gene=gene.names[i],pval=pval)
    if(pval<pval.cut){
      sub.result$DE='yes' } else {sub.result$DE='no'
      }
    results=rbind(results,sub.result)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  de.gene<-as.vector(results[results[,2]<pval.cut,1])
  de.ts<-data.exp[,which(colnames(data.exp) %in% de.gene)]
  return(list(results=results,de.gene=de.gene,de.ts=de.ts))
}
