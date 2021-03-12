SFNG <- function(Nodes, mlinks, begin) {
  pos <- length(begin)
  
  Net <- matrix(0, Nodes, Nodes)
  Net[1:pos, 1:pos] <- begin
  sumlinks <- sum(Net)
  
  while (pos < Nodes) {
    pos <- pos + 1
    linkage <- 0
    while (linkage != mlinks) {
      rnode <- ceiling(runif(1) * pos)
      deg <- sum(Net[, rnode]) * 2
      rlink <- runif(1) * 1
      
      if (rlink < (deg / sumlinks) &&
          Net[pos, rnode] != 1 && Net[rnode, pos] != 1) {
        Net[pos, rnode] <- 1
        Net[rnode, pos] <- 1
        linkage <- linkage + 1
        sumlinks <- sumlinks + 2
      }
    }
  }
  
  # disorder genes
  SFNet <- Net
  ID <- sample(1:Nodes , Nodes)
  SFNet <- SFNet[ID, ID]
  return(SFNet)
}
