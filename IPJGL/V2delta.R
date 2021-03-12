V2delta <- function(theta,V) {
  d <- solve(sqrt(diag(diag(theta))))
  delta <- -d %*% (V+t(V)) %*% d
  diag(delta)<-0
  return(delta)
}