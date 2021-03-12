keep.largest.N <- function(X,N){
  N <- 2*N
  temp<-order(abs(X),decreasing = TRUE)[1:N]
  th<-abs(X[temp[N]])
  X[abs(X)<th]<-0
  return(X)

}
