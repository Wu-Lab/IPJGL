rm(list = ls(all = TRUE))
source('./source.scripts.R')

data <- generate.data(
  n1,
  n2,
  p,
  rates.connect,
  rates.drop,
  m.pert,
  umin = 0.3,
  umax = 0.8,
  diffmode
)
delta.true <- data$delta.true

##### IPJGL #####
results.IPJGL <-
  IPJGL(
    data = data
  )
delta.IPJGL <- theta2partial(results.IPJGL$Z1)-theta2partial(results.IPJGL$Z2)

keep.edges <- p
p1 <- showmatrix(keep.largest.N(delta.true,keep.edges),main = 'truth')
p2 <- showmatrix(keep.largest.N(delta.IPJGL,keep.edges),main = 'IPJGL')
