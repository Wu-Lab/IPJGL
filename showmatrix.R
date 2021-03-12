showmatrix <- function(theta, thresh = 1e-5, main) {
  library('pheatmap')
  theta <- theta - diag(diag(theta))
  range <- seq(-max(abs(theta)), max(abs(theta)), length.out = 100)
  colorbar = colorRampPalette(c("blue", "white", "red"))(100)
  p <- pheatmap(
    theta,
    color = colorbar,
    breaks = range,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = main
  )
  return(p)
}
