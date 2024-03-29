
# cell_abundance_hypergeo_test.R 
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(purrr)
  library(ggplot2)
})

test_for_enrichment <- function(Sample, Cluster, q, m, n, k) {
  phyper(q, m, n, k, lower.tail = FALSE)
}

cell_abundance_hypergeo_test <- function(scrna, compartment) {
  df <- table(scrna$orig.ident, scrna$seurat_clusters) %>% as.data.frame()
  colnames(df) <- c("Sample", "Cluster", "Count")
  df$Sample <- factor(df$Sample, levels = c("Sham", "ICI", "BT", "BTICI",
                                            "2ICI", "8ICI", "20ICI"))
  total_cells <- group_by(df, Sample) %>% summarize(m = sum(Count))
  df <- arrange(df, Cluster, Sample) %>% 
          tidyr::pivot_wider(names_from = Cluster, values_from = Count)
  counts <- as.matrix(df[,-1])
  rownames(counts) <- df$Sample
  
  #q = # of sample-specific cells in cluster
  #m = total # of cells from sample
  #n = total all cells - m
  #k = # of cells in cluster
  # build df that contains that cell counts per sample per cluster
  hyper_df <- table(scrna$orig.ident, scrna$seurat_clusters) %>% as.data.frame()
  colnames(hyper_df) <- c("Sample", "Cluster", "q")
  hyper_df <- left_join(hyper_df, total_cells, by = "Sample")
  hyper_df$n <- sum(total_cells$m) - hyper_df$m
  cluster_totals <- select(hyper_df, Cluster, q) %>% group_by(Cluster) %>%
                      summarise(k = sum(q))
  hyper_df <- left_join(hyper_df, cluster_totals, by = "Cluster")
  # to get P[X>= x]
  hyper_df$q <- hyper_df$q - 1
  
  hyper_df$`P(X>=x)` <- unlist(purrr::pmap(hyper_df, test_for_enrichment))
  hyper_df$fdr <- p.adjust(hyper_df$`P(X>=x)`, method = "BH")
  
  fid <- paste0(compartment, ".cell_abundance_hypergeo_test.results.xlsx")
  openxlsx::write.xlsx(hyper_df, file = fid)

  # Visualize hypergeo p-val results as a heatmap
  df <- select(hyper_df, Sample, Cluster, fdr)
  df$log10p <- -log10(df$fdr)
  
  # define heatmap color scale
  mat_scale <- seq(floor(min(df$log10p, na.rm = TRUE)), 
                   ceiling(max(df$log10p, na.rm = TRUE)), 
                   length.out = 20)
  col_fun <- circlize::colorRamp2(mat_scale, 
                                  pals::brewer.ylorrd(n = length(mat_scale)),
                                  space = "LAB")
  
  df <- select(df, -fdr)
  df <- tidyr::pivot_wider(df, names_from = "Sample", values_from = "log10p")
  mat <- as.matrix(select(df, Sham, ICI, BT, BTICI, `2ICI`, `8ICI`, `20ICI`))
  rownames(mat) <- df$Cluster
  colnames(mat) <- c("Sham", "ICI", "BT", "BT+ICI", "2 Gy + ICI", 
                     "8 Gy + ICI", "20 Gy + ICI")
  
  ht <- Heatmap(mat, name = "-Log10(FDR)", col = col_fun, 
                cluster_rows = FALSE, cluster_columns = FALSE, 
                show_row_names = TRUE, row_names_side = "left",
                row_names_gp = gpar(fontsize = 6), 
                row_title_gp = gpar(fontsize = 8),
                row_title = paste(compartment, "Cluster"), 
                show_column_names = TRUE, column_names_side = "bottom", 
                column_names_gp = gpar(fontsize = 6), 
                column_title_gp = gpar(fontsize = 8),
                #column_title = "Condition", 
                cell_fun = function(j,i,x,y,w,h,col) {
                  if (mat[i,j] > 10) {
                    grid.text("**", x,y, hjust = 0.5, vjust = 0.5, 
                              gp = gpar(fontsize = 6))
                  } else if(mat[i,j] > 5) {
                    grid.text("*", x,y, hjust = 0.5, vjust = 0.5, 
                              gp = gpar(fontsize = 6))
                  }
                }, 
                height = unit(6, "cm"), width = unit(5.25, "cm"), 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 7), 
                                            labels_gp = gpar(fontsize = 6), 
                                            grid_height = unit(1, "cm"), 
                                            grid_width = unit(0.25, "cm")))
  return(ht)
}
