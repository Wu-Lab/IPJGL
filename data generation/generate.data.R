generate.data <-
  function(n1,
           n2,
           p,
           rate.connect = 0.1,
           rate.drop = 0.8,
           m.pert = 10,
           umin = 0.3,
           umax = 0.8,
           diffmode) {
    source('./data generation/SFNG.R')
    source('./data generation/theta2partial.R')
    library(MASS)


    net.structure <- SFNG(p, 2, 1)
    diag(net.structure) <- 0



    dense = matrix(runif(p * p), p)
    dense = (dense - 0.5) / 0.5 * (umax - umin) + umin * sign(dense - 0.5)
    dense[upper.tri(dense)] <- 0
    dense <- dense + t(dense)

    theta1 <- net.structure * dense
    theta2 <- theta1

    all.index <- 1:p

    #group genes
    score <- rep(0, p)
    for (i in 1:p) {
      score[i] <- sum(net.structure[, i])
    }
    index.sort <- order(score, decreasing = TRUE)

    index.vital <- index.sort[1:ceiling(0.2 * p)]
    index.moderate <-
      index.sort[(ceiling(0.2 * p) + 1):ceiling(0.4 * p)]
    index.original <- index.sort[(ceiling(0.4 * p) + 1):p]

    # select diffgene
    if (diffmode == 1) {
      index.pert <- sample(index.vital, m.pert, replace = FALSE)
    } else if (diffmode == 2) {
      index.pert <- sample(index.moderate, m.pert, replace = FALSE)

    } else if (diffmode == 3) {
      index.pert <- sample(index.original, m.pert, replace = FALSE)
    } else if (diffmode == 4) {
      index.pert <- sample(index.vital, m.pert / 2, replace = FALSE)
      index.pert <-
        c(index.pert,
          sample(index.original, m.pert / 2, replace = FALSE))
    } else if (diffmode == 5) {
      index.pert <- sample(index.vital, m.pert / 10, replace = FALSE)
      index.pert <-
        c(index.pert,
          sample(index.moderate, m.pert / 10 * 4, replace = FALSE))
      index.pert <-
        c(index.pert,
          sample(index.original, m.pert / 10 * 5, replace = FALSE))
    }


    #  pertubation
    for (pert in index.pert) {
      connect.now <- which(net.structure[pert, ] != 0)
      disconnect.now <- all.index[-c(connect.now, pert)]

      disconnect.index <- which(runif(p) < rate.drop)
      disconnect.index <- intersect(connect.now, disconnect.index)

      theta2[pert, disconnect.index] <- 0
      theta2[disconnect.index, pert] <- 0

      connect.index <- which(runif(p) < rate.connect)
      connect.index <- intersect(disconnect.now, connect.index)

      l3 <- length(connect.index)
      connect.dense <- runif(l3)
      connect.dense <-
        (connect.dense - 0.5) / 0.5 * umin + (umax - umin) * sign(connect.dense -
                                                                    0.5)
      theta2[pert, connect.index] <- connect.dense
      theta2[connect.index, pert] <- connect.dense
    }
    if (sum((theta1 - theta2) != 0) == 0) {
      for (pert in index.pert) {
        connect.now <- which(net.structure[pert, ] != 0)
        disconnect.now <- all.index[-c(connect.now, pert)]

        disconnect.index <- which(runif(p) < rate.drop)
        disconnect.index <-
          intersect(connect.now, disconnect.index)

        theta2[pert, disconnect.index] <- 0
        theta2[disconnect.index, pert] <- 0

        connect.index <- which(runif(p) < rate.connect)
        connect.index <- intersect(disconnect.now, connect.index)

        l3 <- length(connect.index)
        connect.dense <- runif(l3)
        connect.dense <-
          (connect.dense - 0.5) / 0.5 * umin + (umax - umin) * sign(connect.dense -
                                                                      0.5)
        theta2[pert, connect.index] <- connect.dense
        theta2[connect.index, pert] <- connect.dense
      }
    }

    eig1 <- eigen(theta1)$values
    eig2 <- eigen(theta2)$values
    eig.min <- min(c(eig1, eig2))
    eig.diag <- (abs(eig.min) + 0.1) * diag(p)
    theta1 <- theta1 + eig.diag
    theta2 <- theta2 + eig.diag

    X1 <- mvrnorm(n1, rep(0, p), solve(theta1))
    X2 <- mvrnorm(n2, rep(0, p), solve(theta2))

    sigma1 <- cov(X1)
    sigma2 <- cov(X2)

    thetap <- lapply(list(theta1, theta2), theta2partial)
    delta.true <- thetap[[1]] - thetap[[2]]

    return(
      list(
        delta.true = delta.true,
        theta1 = theta1,
        theta2 = theta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        X1 = X1,
        X2 = X2,
        index.pert = index.pert
      )
    )
  }
