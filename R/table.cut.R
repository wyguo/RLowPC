#' Cut the confusion table at threshods
#'
#' The confusion table output from \code{\link{table.evaluate}} is cut at a threshold of top weighted edges, FPR or TPR. The \code{cutoff} in the function must match
#' to the \code{cutat} metrices.
#' @param input.table the confusion table output from \link{table.evalution}.
#' @param cutoff a numeric value of threshold for "top", "fpr" or "tpr".
#' @param cutat the options are "top", "fpr" and "tpr", which indicate to cut the confusion table at top weighted edges, FPR or TPR.
#' @return \code{table.cut} returns a list of results cantains the table after applying cutoff, the pAUROC value, pAUPR value, precesion for the remained edges, the corresponding top, fpr, tpr cutoffs
#' and the statistical measures after applying function \link{confusion}.
#' @export
table.cut<-function(input.table,cutoff=1000,cutat='top'){
  if(cutat=='top'){
    input.table$idx=input.table$tp+input.table$fp
    table.output<-input.table[which(input.table$idx<=cutoff),1:5]
  }
  if(cutat=='tpr'){
    if(cutoff>1)
      stop('tpr must less than 1')
    input.table$idx<-input.table$tp/(input.table$tp+input.table$fn)
    table.output<-input.table[input.table$idx<=cutoff,1:5]
  }
  if(cutat=='fpr'){
    if(cutoff>1)
      stop('fpr must less than 1')
    input.table$idx<-input.table$fp/(input.table$fp+input.table$tn)
    table.output<-input.table[input.table$idx<=cutoff,1:5]
  }
  roc.output<-minet::auc.roc(table.output)
  pr.output<-minet::auc.pr(table.output)
  top.output<-table.output$tp+table.output$fp

  #######
  table.conf<-confusion(table.output)
  table.conf$id<-table.conf$tp+table.conf$fp
  #table.conf<-table.conf[!duplicated(table.conf$id),]
  precision.output<-table.conf$precision
  fpr.output<-table.conf$fpr
  tpr.output<-table.conf$tpr
  results=list(table.output=table.output,
               roc=roc.output,pr=pr.output,precision=precision.output,
               fpr=fpr.output,tpr=tpr.output,top=top.output,
               table.conf=table.conf)
  return(results)
}
