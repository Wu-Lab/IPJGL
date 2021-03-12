soft.tauV <- function(X, theta1, theta2, lambda2, rho, normtype) {
  p <- length(X[, 1])
  V <- matrix(0, p, p)
  w <- matrix(0, 1, p)
  for (i in 1:p) {
    w[i] <- norm(theta1[, i], '2') ^ 2 + norm(theta2[, i], '2') ^ 2
  }
  
  if (normtype == '1') {
    for (i in 1:p) {
      V[, i] <-
        sign(X[, i]) * (pmax(0, abs(X[, i]) - lambda2 * w[i])) / (2 *
                                                                    rho)
      
    }
  } else if (normtype == '2') {
    for (i in 1:p) {
      if (norm(X[, i], '2') > lambda2 * w[i]) {
        V[, i] <-
          (norm(X[, i], '2') - lambda2 * w[i]) / (2 * rho * norm(X[, i], '2')) *
          X[, i]
      }
    }
  }
  return(V)
}
