soft.tau1 <- function(X, lambda) {
  Y <- sign(X) * (pmax(abs(X) - lambda, 0))
  return(Y)
}
