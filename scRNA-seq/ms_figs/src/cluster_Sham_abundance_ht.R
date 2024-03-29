
#cluster_Sham_abundance_ht.R

# Load libraries silently
suppressPackageStartupMessages({
  library(ComplexHeatmap)
})

cluster_Sham_abundance_ht <- function(scrna, compartment) {
  sample_counts <- table(scrna$orig.ident)
  
  df <- table(scrna$orig.ident, scrna$seurat_clusters) %>% as.data.frame()
  colnames(df) <- c("Sample", "Cluster", "Count")
  df <- mutate(df, Norm_count = Count/sample_counts[Sample]) %>% arrange(Cluster)
  sham_norm_count <- filter(df, Sample == "Sham") %>% pull(Norm_count)
  names(sham_norm_count) <- levels(df$Cluster)
  
  df <- mutate(df, log2FC = log2(Norm_count/sham_norm_count[Cluster]))
  df$Sample <- factor(df$Sample, levels = c("Sham", "ICI", "BT", "BTICI", 
                                            "2ICI", "8ICI", "20ICI"))
  df <- arrange(df, Cluster, Sample)
  mat <- matrix(df$log2FC, ncol = length(levels(df$Sample)), byrow = TRUE)
  colnames(mat) <- levels(df$Sample)
  colnames(mat) <- c("Sham", "ICI", "BT", "BT+ICI", "2 Gy + ICI", 
                     "8 Gy + ICI", "20 Gy + ICI")
  rownames(mat) <- levels(df$Cluster)
  h <- nrow(mat)/2
  # remove sham form mat
  mat <- mat[,-1]
  
  mat_scale <- seq(floor(min(mat)), ceiling(max(mat)), by = 2)
  col_fun <- circlize::colorRamp2(mat_scale, space = "LAB",
                                  pals::coolwarm(n = length(mat_scale)))
  ht <- Heatmap(mat, name = "log2\n(Sample/Sham)", col = col_fun, 
                cluster_rows = TRUE, cluster_columns = FALSE, 
                show_row_names = TRUE, row_names_side = "left", 
                show_column_names = TRUE, column_names_side = "bottom", 
                row_title = paste(compartment, "Cluster"), 
                row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                column_title = "Cluster Abundance Relative to Sham", 
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10), 
                row_names_gp = gpar(fontsize = 10),
                width = unit(5.5, "cm"), height = unit(h, "cm"), 
                border = TRUE, border_gp = gpar(col = "grey40", lwd = 0.5), 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
                                            labels_gp = gpar(fontsize = 8)))
  
  return(ht)
}