#' Statistical derivations of confusinon table
#'
#' Calculate statistical measures of the performance of binary classification test from the output confusion matrix \code{\link{table.evaluate}}
#'
#' @param input.table the output confusion table from \code{\link{table.evaluate}}
#' @details
#' true positive: \eqn{tp}; false positive: \eqn{fp}; true negative: \eqn{tn}; false negative: \eqn{fn};\cr
#' positives in reference network: \eqn{p}; negatives in reference network: \eqn{n};\cr
#' true positive rate: \eqn{tpr=recall=\frac{tp}{tp+fn}}; false positive rate: \eqn{fpr=\frac{fp}{fp+tn}}; \cr
#' true negative rate: \eqn{tnr=\frac{tn}{tn+fp}}; false negative rate: \eqn{fpr=\frac{fn}{fn+tp}};\cr
#' precision: \eqn{precision=\frac{tp}{tp+fp}}; negative predictive value: \eqn{npv=\frac{tn}{tn+fn}};\cr
#' false discovery rate: \eqn{fdr=\frac{fp}{fp+tp}}; accuracy: \eqn{accuracy=\frac{tp+tn}{p+n}};\cr
#' f1 scaore: \eqn{f1=\frac{2tp}{2tp+fp+fn}}; \cr
#' Matthews correlation coefficient: \deqn{mcc=\frac{tp\times tn-fp\times fn}{\sqrt{(tp+fp)\times (tp+fn)\times (tn+fp)\times (tn+fn)}}}
#'
#' @return a table of of statistical measures of performance, see Details.
#' @references Powers DMW: Evaluation: From Precision, Recall and F-Factor to ROC, Informedness, Markedness & Correlation. In. Adelaide, Australia; 2007.
#' @export



confusion<-function(input.table){
  table<-input.table
  tp<-table$tp
  fp<-table$fp
  tn<-table$tn
  fn<-table$fn
  p<-tp+fn
  n<-fp+tn
  tpr<-tp/(tp+fn)
  fpr<-fp/(fp+tn)
  tnr<-tn/(fp+tn)#or specificity(spc)
  fnr<-fn/(fn+tp)#or miss rate
  precision<-tp/(tp+fp)
  recall<-tp/(tp+fn)# or tpr
  npv<-tn/(tn+fn)#negative predictive value
  fdr<-fp/(fp+tp)#false discovry rate
  accuracy<-(tp+tn)/(p+n)
  f1<-2*tp/(2*tp+fp+fn)
  mcc<-(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  result<-data.frame(thrsh=table$thrsh,tp=tp,fp=fp,fn=fn,tn=tn,tpr=tpr,fpr=fpr,tnr=tnr,fnr=fnr,precision=precision,recall=recall,npv=npv,
                     fdr=fdr,accuracy=accuracy,f1=f1,mcc=mcc)
  return(result)
}
