#! Rscript

# Fig5G.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(glmGamPoi)
  library(DESeq2)
  library(Seurat)
  library(ComplexHeatmap)
})

lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

results <- tibble(Cluster = as.character(unique(lymphoid_scrna$seurat_clusters)))
results <- mutate(results, DESeq_LRT = purrr::map(Cluster, 
                                                  function(c) {synergistic_DE(scrna = lymphoid_scrna, 
                                                                              cluster = c)}))

# compile results into single df for each cluster for writing to Excel
results_xlsx <- mutate(results, compiled = purrr::map(DESeq_LRT, annotate_EBRT_padj))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Column info")
df <- data.frame(Column = colnames(results_xlsx$compiled[[1]]), 
                 Description = c("MGI gene symbol", "average expression across Sham, ICI, BT, BTICI", 
                                 "of BTICI relative to Sham", "ignore", "ignore", 
                                 "FDR of BTICI interaction effect", 
                                 "% of Sham cells with reads aligning to the gene", 
                                 "% of ICI cells with reads aligning to the gene", 
                                 "% of BT cells with reads aligning to the gene", 
                                 "% of BTICI cells with reads aligning to the gene", 
                                 "FDR of 2ICI interaction effect", 
                                 "% of 2ICI cells with reads aligning to the gene",
                                 "FDR of 8ICI interaction effect", 
                                 "% of 8ICI cells with reads aligning to the gene", 
                                 "FDR of 20ICI interaction effect", 
                                 "% of 20ICI cells with reads aligning to the gene"))
openxlsx::writeData(wb, sheet = "Column info", x = df)

for (n in 1:nrow(results_xlsx)) {
  cluster = results_xlsx$Cluster[n]
  openxlsx::addWorksheet(wb, sheetName = paste("Cluster", cluster))
  df <- results_xlsx$compiled[[n]]
  openxlsx::writeData(wb, sheet = paste("Cluster", cluster), 
                      x = df, keepNA = TRUE)
}

result_stats <- tibble(Cluster = results_xlsx$Cluster, Compartment = "lymphoid")
result_stats$data <- purrr::map(results_xlsx$compiled, count_hits)
result_stats <- tidyr::unnest(result_stats, cols = data)
openxlsx::addWorksheet(wb, sheetName = "Results Statistics")
openxlsx::writeData(wb, sheet = "Results Statistics", 
                    x = result_stats, keepNA = TRUE)

fid <- paste0(fid.prefix, ".lymphoid_clusters.xlsx")
openxlsx::saveWorkbook(wb, file = fid, overwrite = TRUE)

# visualize cluster 11 results in heatmap
# selecting cluster from original results obj.
results <- filter(results, Cluster == "11")

purrr::walk2(results$Cluster, results$DESeq_LRT, 
             plot_interaction_hits, 
             scrna = lymphoid_scrna)
ht <- plot_interaction_hits("11", results$DESeq_LRT[[1]], scrna = )
fid <- "Fig5G.lymphoid.cluster_11.interaction_genes.heatmap.png"
png(fid, height = 27, width = 11.5, res = 200, units = "cm")
draw(ht)
invisible(dev.off())

