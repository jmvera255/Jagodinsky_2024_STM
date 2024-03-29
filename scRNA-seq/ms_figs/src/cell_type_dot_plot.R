
# cell_type_dot_plot.R
cell_type_dot_plot <- function(scrna){
  # dotplot of cell types per cluster
  df <- as.data.frame(table(scrna$Alt_cell_type, scrna$seurat_clusters))
  cluster_sums <- as.data.frame(table(scrna$seurat_clusters))
  df <- left_join(df, cluster_sums, by = c("Var2" = "Var1"))
  colnames(df) <- c("cell_type", "cluster", "cell_count", "cluster_total_cells")
  df$pct_type <- df$cell_count/df$cluster_total_cells
  
  p <- ggplot(df, aes(y = cluster, x = cell_type, 
                      size = pct_type*100, color = pct_type*100)) +
    geom_point() +
    scale_color_gradient(high = "mediumpurple", low = "mistyrose2") +
    ggtitle("Cell Type Abundance by Cluster") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  p$labels$size <- "Fraction of cluster"
  
  return(p)
}

