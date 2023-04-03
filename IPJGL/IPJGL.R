IPJGL <- function(data,
                  lambda1 = 4.6,
                  lambda2 = 5.2,
                  err.threshold = 1e-4,
                  step.max = 1e2,
                  truncate = 1e-5,
                  normtype = '2') {
  source('./IPJGL/soft.tau1.R')
  source('./IPJGL/soft.tauV.R')
  source('./IPJGL/solve.theta.R')
  source('./IPJGL/Dnorm.R')
  library('progress')
  # calculate n1,n2,p
  sigma1 <- cov(data$X1)
  sigma2 <- cov(data$X2)
  n1 <- nrow(data$X1)
  n2 <- nrow(data$X2)
  p <- ncol(data$X1)

  I <- diag(p)

  theta1 <- I

  theta2 <- I

  theta1.old <- I

  theta2.old <- I

  Z1 <- I

  Z2 <- I

  V <- I

  W <- I

  O <- matrix(0, p, p)

  F <- O

  G <- O

  Q1 <- O

  Q2 <- O




  rho0 <- 50

  miu <- 5

  rho.max <- 1e10

  step.max <- 1e2

  rho <- rho0


  rho.steps <- floor(log10(rho.max) / log10(miu))

  pb <- progress_bar$new(total = rho.steps * step.max)

  ##### ADMM #####
  for (j in 1:rho.steps) {
    for (step in 1:step.max) {
      pb$tick()

      theta1 <-
        solve.theta((n1 * sigma1 + F + Q1) - rho * (theta2 + (V + W) + Z1),
                    V,
                    rho,
                    n1,
                    lambda2,
                    normtype)

      theta2 <-
        solve.theta((n2 * sigma2 - F + Q2) - rho * (theta1 - (V + W) + Z2),
                    V,
                    rho,
                    n2,
                    lambda2,
                    normtype)

      V <-
        soft.tauV(F - G + rho * (theta1 - theta2 - W + t(W)),
                  theta1,
                  theta2,
                  lambda2,
                  rho,
                  normtype)

      Z1 <- soft.tau1(theta1 + Q1 / rho, lambda1 / rho)

      Z2 <- soft.tau1(theta2 + Q2 / rho, lambda1 / rho)


      W <-
        0.5 * (t(V) - V + theta1 - theta2) + 0.5 / rho * (F + t(G))
      F <- F + rho * (theta1 - theta2 - (V + W))
      G <- G + rho * (V - t(W))

      Q1 <- Q1 + rho * (theta1 - Z1)

      Q2 <- Q2 + rho * (theta2 - Z2)

      err <-
        max(
          norm(as.matrix(theta1 - theta1.old), 'f') / norm(theta1.old, 'f'),
          norm(as.matrix(theta2 - theta2.old), 'f') / norm(theta2.old, 'f')
        )

      if (err <= err.threshold) {
        break
      }

      theta1.old <- theta1

      theta2.old <- theta2

    }
    rho <- rho * miu

  }
  cat(paste0('\n algorithm is done! relative error=', round(err, 8)))

  theta1 <- (theta1 + t(theta1)) / 2
  theta2 <- (theta2 + t(theta2)) / 2
  Z1 <- (Z1 + t(Z1)) / 2
  Z2 <- (Z2 + t(Z2)) / 2
  theta1[abs(theta1) < truncate] <- 0
  theta2[abs(theta2) < truncate] <- 0
  Z1[abs(Z1) < truncate] <- 0
  Z2[abs(Z2) < truncate] <- 0
  V[abs(V) < truncate] <- 0

  return(
    list(
      theta1 = theta1,
      theta2 = theta2,
      Z1 = Z1,
      Z2 = Z2,
      V = V
    )
  )
}
