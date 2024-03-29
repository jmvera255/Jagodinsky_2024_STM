#! Rscript

# 06_scLANE.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(scran)
  library(purrr)
  library(scater)
  library(scLANE)
  library(ggplot2)
  library(scaffold)
  library(SingleCellExperiment)
  library(Seurat)
})

source("src/functions_simulation.R")

select <- dplyr::select
filter <- dplyr::filter

spe <- qs::qread("BTD3-C1.seurat_spatial.qs")
spe <- DietSeurat(spe, assay = "SCT")

# Get BT distance ccoldata
fid <- "BTD3-C1.BT_dose_annotations.qs"
BT_coldata <- qs::qread(fid)

# filter to spots in seurat obj
BT_coldata <- filter(BT_coldata, barcode %in% colnames(spe))
rownames(BT_coldata) <- BT_coldata$barcode

# add BT distance to spe
spe[['dist_um']] <- BT_coldata[colnames(spe), "dist_um"]

dist_max <- max(spe$dist_um)
spe[['dist_norm']] <- spe$dist_um/dist_max
order_df <- data.frame(X = spe$dist_norm)
dist_df <- data.frame(X = spe$dist_um)

# there shouldn't be any offset since SCT has normalized for seq depth
cell_offset <- rep(1, ncol(spe))
names(cell_offset) <- colnames(spe)

de_test_glm <- testDynamic(spe, 
                           pt = order_df, 
                           genes = rownames(spe),
                           size.factor.offset = cell_offset, 
                           n.cores = 6, 
                           track.time = TRUE)

fid <- "BTD3-C1.scLANE.test_glm.qs"
qs::qsave(de_test_glm, fid)

# apply statistical test to find dynamically expressed genes
de_res_glm <- getResultsDE(de_test_glm) %>% 
  select(Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall)

fid <- "BTD3-C1.scLANE.res_glm.qs"
qs::qsave(de_res_glm, fid)

fid <- "BTD3-C1.scLANE.res_glm.results.xlsx"
openxlsx::write.xlsx(x = de_res_glm, file = fid, overwrite = TRUE)


# END