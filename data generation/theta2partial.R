theta2partial <- function(theta) {
  d <- solve(sqrt(diag(diag(theta))))
  partial <- -d %*% theta %*% d
  diag(partial)<-1
  return(partial)
}