Dnorm <- function(V, normtype)
{
  p = length(V[1, ])
  
  DqV = matrix(0, p, p)
  
  for (i in 1:p)
  {
    DqV[i, i] = norm(as.matrix(V[ , i]), normtype)
  }
  return(DqV)
}