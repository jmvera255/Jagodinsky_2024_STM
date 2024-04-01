#! Rscript

# Fig4F.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(glmGamPoi)
  library(DESeq2)
  library(Seurat)
  library(ComplexHeatmap)
})

source("src/Fig4F.functions.R")

# load clustered myeloid compartment object
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")

# subset for neutrophils cluster
myeloid_integrated <- subset(myeloid_integrated, subset = seurat_clusters == "8")


results <- tibble(Cluster = as.character(unique(myeloid_integrated$seurat_clusters)))
results <- mutate(results, DESeq_LRT = purrr::map(Cluster, 
                                                  function(c) {condition_DE(scrna = myeloid_integrated,
                                                                            cluster = c)}))

# compile results into single df for each cluster for writing to Excel
results <- mutate(results, compiled = purrr::map(DESeq_LRT, combine_cluster_results))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Column info")
df <- data.frame(Column = colnames(results$compiled[[1]]), 
                 Description = c("MGI gene symbol", 
                                 "% of Sham cells with reads aligning to the gene",
                                 "log2 fold change ICI relative to Sham", 
                                 "raw pvalue of likelihood ratio test", "FDR of likelihood ratio test", 
                                 "% of ICI cells with reads aligning to the gene", 
                                 "log2 fold change BT relative to Sham", 
                                 "raw pvalue of likelihood ratio test", "FDR of likelihood ratio test", 
                                 "% of BT cells with reads aligning to the gene", 
                                 "log2 fold change BTICI relative to Sham", 
                                 "raw pvalue of likelihood ratio test", "FDR of likelihood ratio test", 
                                 "% of BTICI cells with reads aligning to the gene",
                                 "log2 fold change 2ICI relative to Sham", 
                                 "raw pvalue of likelihood ratio test", "FDR of likelihood ratio test", 
                                 "% of 2ICI cells with reads aligning to the gene",
                                 "log2 fold change 8ICI relative to Sham", 
                                 "raw pvalue of likelihood ratio test", "FDR of likelihood ratio test", 
                                 "% of 8ICI cells with reads aligning to the gene", 
                                 "log2 fold change 20ICI relative to Sham", 
                                 "raw pvalue of likelihood ratio test", "FDR of likelihood ratio test", 
                                 "% of 20ICI cells with reads aligning to the gene"))
openxlsx::writeData(wb, sheet = "Column info", x = df)
for (n in 1:nrow(results)) {
  cluster = results$Cluster[n]
  openxlsx::addWorksheet(wb, sheetName = paste("Cluster", cluster))
  df <- results$compiled[[n]]
  openxlsx::writeData(wb, sheet = paste("Cluster", cluster), 
                      x = df, keepNA = TRUE)
}

result_stats <- tibble(Cluster = results$Cluster, Compartment = "myeloid")
result_stats$data <- purrr::map(results$compiled, count_hits)
result_stats <- tidyr::unnest(result_stats, cols = data)
openxlsx::addWorksheet(wb, sheetName = "Results Statistics")
openxlsx::writeData(wb, sheet = "Results Statistics", 
                    x = result_stats, keepNA = TRUE)

fid <- "Fig4F.DESeq_results.xlsx"
openxlsx::saveWorkbook(wb, file = fid, overwrite = TRUE)

### Visualize pro-/anti-tumor gene expression heatmap in neutrophils
antitumor <- c("Icam1", "Cd177", "Itgam", "Tnf", "C3")
protumor <- c("Cd274", "Cd14", "Cd101", "Sirpa", "Siglecf", "Ccl17", "Adgrg5")

results <- filter(results, Cluster == "8")
df <- results$compiled[[1]]

antitumor <- antitumor[which(antitumor %in% df$Gene)]
protumor <- protumor[which(protumor %in% df$Gene)]

df <- filter(df, Gene %in% c(antitumor, protumor))

lfc_mat <- as.matrix(select(df, matches("LFC_")))
rownames(lfc_mat) <- df$Gene
lfc_mat <- lfc_mat[c(protumor,antitumor),]
colnames(lfc_mat) <- stringr::str_remove(colnames(lfc_mat), "LFC_")
colnames(lfc_mat) <- stringr::str_replace(colnames(lfc_mat), 
                                          "BTICI", "BT+ICI")
colnames(lfc_mat) <- stringr::str_replace(colnames(lfc_mat), 
                                          "2ICI", "2 Gy + ICI")
colnames(lfc_mat) <- stringr::str_replace(colnames(lfc_mat), 
                                          "8ICI", "8 Gy + ICI")
colnames(lfc_mat) <- stringr::str_replace(colnames(lfc_mat), 
                                          "20ICI", "20 Gy + ICI")

fdr_mat <- as.matrix(select(df, matches("padj_")))
colnames(fdr_mat) <- stringr::str_remove(colnames(fdr_mat), "padj_")
colnames(fdr_mat) <- stringr::str_replace(colnames(fdr_mat), 
                                          "BTICI", "BT+ICI")
colnames(fdr_mat) <- stringr::str_replace(colnames(fdr_mat), 
                                          "2ICI", "2 Gy + ICI")
colnames(fdr_mat) <- stringr::str_replace(colnames(fdr_mat), 
                                          "8ICI", "8 Gy + ICI")
colnames(fdr_mat) <- stringr::str_replace(colnames(fdr_mat), 
                                          "20ICI", "20 Gy + ICI")
rownames(fdr_mat) <- df$Gene
# make sure fdr and lfc mats are in same order
fdr_mat <- fdr_mat[rownames(lfc_mat),colnames(lfc_mat)]
fdr_mat[is.na(fdr_mat)] <- 1

x <- max(abs(c(lfc_mat)), na.rm = TRUE) %>% ceiling()
mat_scale <- seq(-1*x, x, by = 1)
col_fun <- circlize::colorRamp2(mat_scale, 
                                pals::coolwarm(n = length(mat_scale)),
                                space = "LAB")
rowsplit <- c(rep("Protumor", times = length(protumor)), 
              rep("Antitumor", times = length(antitumor)))
rowsplit <- factor(rowsplit, levels = c("Protumor", "Antitumor"))

ht <- Heatmap(lfc_mat, col = col_fun, name = "log2\nfold\nchange", 
              row_split = rowsplit, row_title = levels(rowsplit), 
              cluster_rows = FALSE, 
              row_names_gp = gpar(fontsize = 7), 
              row_title_gp = gpar(fontsize = 9),
              show_column_names = TRUE, column_names_side = "bottom",
              cluster_columns = TRUE, column_dend_side = "top",
              column_dend_height = unit(0.35, "cm"),
              column_title = "Neutrophils (Myeloid Cluster 8)",
              column_names_gp = gpar(fontsize = 7), 
              column_title_gp = gpar(fontsize = 9),
              cell_fun = function(j,i,x,y,w,h,col) {
                if (fdr_mat[i,j] < 0.0001) {
                  grid.text("**", x,y, hjust = 0.5, vjust = 0.5, 
                            gp = gpar(fontsize = 6))
                } else if(fdr_mat[i,j] < 0.05) {
                  grid.text("*", x,y, hjust = 0.5, vjust = 0.5, 
                            gp = gpar(fontsize = 6))
                }
              }, 
              height = unit(4.5, "cm"), width = unit(5.25, "cm"),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 6), 
                                          labels_gp = gpar(fontsize = 5), 
                                          grid_height = unit(0.85, "cm"), 
                                          grid_width = unit(0.25, "cm")))

fid <- "Fig4F.neutrophil.gene_expression.heatmap.png"
png(fid, width = 8.75, height = 7, units = "cm", res = 300)
draw(ht)
invisible(dev.off())


# END