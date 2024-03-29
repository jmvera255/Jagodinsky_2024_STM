#! Rscript

# 05_Seurat_transferData.R


# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(SpatialExperiment)
  library(glmGamPoi)
})

source("src/Seurat_transferData.functions.R")
`%notin%` <- Negate(`%in%`)

# read in gene_meta so ensemble ids can be converted to symbols to match ids used in scRNA
gene_meta <- qs::qread("2022-11-21_spotclean.gene_meta.qs")
duplicates <- table(gene_meta$mgi_symbol) %>% as.data.frame() %>% 
  filter(Freq > 1) %>% pull(Var1) %>% as.character()
gene_meta <- filter(gene_meta, mgi_symbol %notin% duplicates)
rownames(gene_meta) <- gene_meta$ensembl_gene_id

# get scrna seurat object list prepped in 04_prep_scRNA_for_TransferData.R
fid <- "scrna_list.TransferData_input.qs"
scrna_list <- qs::qread(fid)
# select BT+ICI scrna-seq data
BTICI_scrna <- scrna_list$BTICI

# get BTD3-C1 clustered spe object
fid <- "BTD3-C1.clusters.spe.qs"
spe <- qs::qread(fid)

# convert BTD3-C1 SpatialExperiment object to Seurat Spatial object
seurat_spe <- convert_to_seurat_spatial(spe)

# cell cylce scoring of spatial data
# get cell cycle marker genes
fid <- "../supporting_data/cell_cycle_genes.qs"
cell_cycle_genes <- qs::qread(fid)

seurat_spe <- CellCycleScoring(seurat_spe, s.features = unlist(cell_cycle_genes$S), 
                        g2m.features = unlist(cell_cycle_genes$`G2/M`))
seurat_spe$CC.Difference <- seurat_spe$S.Score - seurat_spe$G2M.Score  

# SCTransform seurat_spe
N <- length(seurat_spe$orig.ident)
seurat_spe <- SCTransform(seurat_spe, vst.flavor = "v2", 
                          verbose = FALSE, assay = "Spatial",
                          vars.to.regress = c("CC.Difference"), 
                          variable.features.n = 2500, ncells = N)
seurat_spe <- RunPCA(seurat_spe, verbose = FALSE, npcs = 50)

# perform reference mapping
anchors <- FindTransferAnchors(reference = BTICI_scrna, query = seurat_spe, 
                               normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = BTICI_scrna$Alt_cell_type, 
                                  prediction.assay = TRUE, 
                                  weight.reduction = seurat_spe[["pca"]], 
                                  dims = 1:30)
# add cell type predictions as assay to seurat_spe
seurat_spe[["scRNA_cell_type"]] <- predictions.assay

# save new Seurat Spatial object with cell type predictions
fid <- "BTD3-C1.seurat_spatial.qs"
qs::qsave(seurat_spe, fid)


# END