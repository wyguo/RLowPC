orderadjmatrix<-function(mat,order){
  m<-matrix(0,length(order),length(order));
  dimnames(m)<-list(order,order);
  index<-match(colnames(mat),colnames(m));
  m[index,index]<-mat;
  mat<-m;
  return(m)
}