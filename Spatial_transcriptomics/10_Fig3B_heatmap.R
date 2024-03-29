#! Rscript

# 10_Fig3B_heatmap.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  
})

# get Seurat spatial data object with cell type predictions assay made in 05_Seurat_transferData.R
fid <- "BTD3-C1.seurat_spatial.qs"
spe <- qs::qread(fid)

# get BT dose annotation colData
fid <- "BTD3-C1.BT_dose_annotations.qs"
BT_coldata <- qs::qread(fid) %>% select(barcode, BayesSpace_clusters, BT_dose)

# get new CL-based cell type predictions

cl_df <- as.matrix(GetAssayData(spe, 
                                assay = "scRNA_cell_type", 
                                slot = "data"))
cl_df <- as.data.frame(t(cl_df)) %>% select(-max)
# find and remote cell types that are not predicted in the sample
keep <- colSums(cl_df) > 0
cl_df <- cl_df[,keep]
cl_df$barcode <- rownames(cl_df)
cl_df <- left_join(cl_df, BT_coldata, by = "barcode")

# wrangle data into long-format for subsequent collapse into average score
cl_df <- tidyr::pivot_longer(cl_df, cols = matches("^CL"), 
                             names_to = "CL", values_to = "Score")
cl_df <- left_join(cl_df, select(mouse_imm_df_unique, CL, Cell_type), by = "CL")

# immune cell predictions/composition per spot
df <- group_by(cl_df, BT_dose, q2_clusters, Cell_type) %>% 
  summarise(Score = mean(Score), .groups = "drop")
df <- tidyr::pivot_wider(df, names_from = Cell_type, values_from = Score)

# count how many cells per cluster per BT_dose
cell_count <- select(BT_coldata, q2_clusters, BT_dose) %>% 
  group_by(q2_clusters, BT_dose) %>%
  summarise(N = n(), .groups = "drop")

# add counts to df
df <- left_join(df, cell_count, by = c("q2_clusters", "BT_dose"))

# rearrange columns
df <- select(df, q2_clusters, BT_dose, N, everything())

df$q2_clusters <- factor(as.character(df$q2_clusters, levels = sort(unique(as.numeric(df$q2_clusters)))))
df$BT_dose <- factor(df$BT_dose, levels = c("high", "mod", "low"))
df <- arrange(df, q2_clusters, BT_dose)

row_ha <- get_row_ha(select(df, BT_dose, q2_clusters, N))
mat <- as.matrix(select(df, -BT_dose, -q2_clusters, -N))
colnames(mat) <- colnames(df)[-c(1,2,3)]
if (!by_cluster) {
  row_ha@anno_list$Cluster <- NULL
  row_ha@anno_list$Hypergeo <- NULL
  row_ha@anno_size <- row_ha@anno_size[2:3]
  row_ha@width <- unit(13, "mm")
}
plot_predictions_heatmap(mat, row_ha, sample, by_cluster, rsplit = df$q2_clusters)