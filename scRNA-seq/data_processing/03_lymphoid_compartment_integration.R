#! Rscript
# 03_lymphoid_compartment_integration.R
# lymphoid.3000_20_0.8

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# get lymphoid_scrna_list from scrna_qc.R
fid <- "scrna_qc.lymphoid.scrna_list.processed.qs"
lymphoid_scrna_list <- qs::qread(fid)
samples <- names(lymphoid_scrna_list)

# apply standard normalization to data
lymphoid_scrna_list <- purrr::map(lymphoid_scrna_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
  return(x)
})

# perform cell cylce scoring and SCTransform with regression of cell cycle and mito reads
cell_cycle_genes <- qs::qread("../supporting_data/cell_cycle_genes.qs")

lymphoid_scrna_list <- purrr::map(lymphoid_scrna_list, function(scrna) {
  N <- ncol(GetAssayData(scrna, slot = "counts"))
  scrna <- CellCycleScoring(scrna, s.features = unlist(cell_cycle_genes$S), 
                            g2m.features = unlist(cell_cycle_genes$`G2/M`))
  scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score
  scrna <- SCTransform(scrna, vst.flavor = "v2", verbose = FALSE, 
                       vars.to.regress = c("CC.Difference", "percent_mito"), 
                       variable.features.n = 3000, ncells = N)
  scrna <- RunPCA(scrna, verbose = FALSE, npcs = 50)
  return(scrna)
})

# find and prep integration anchors
features <- SelectIntegrationFeatures(lymphoid_scrna_list, nfeatures = 3000)
lymphoid_scrna_list <- PrepSCTIntegration(lymphoid_scrna_list, anchor.features = features)
anchors <- FindIntegrationAnchors(lymphoid_scrna_list, normalization.method = "SCT", 
                                  anchor.features = features)

# integrate samples
scrna_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(scrna_integrated) <- "integrated"

#then perform standard downstream steps
scrna_integrated <- RunPCA(scrna_integrated, verbose = FALSE, npcs = 50)
scrna_integrated <- RunUMAP(scrna_integrated, dims = 1:20)
scrna_integrated <- FindNeighbors(scrna_integrated, dims = 1:20, verbose = FALSE)

# find clusters using Louvain algorithm with multilevel refinement
scrna_integrated <- FindClusters(scrna_integrated, resolution = 0.8, 
                                 random.seed = 1234, algorithm = 2, 
                                 verbose = FALSE)
fid <- "lymphoid_integrated.scrna.qs"
qs::qsave(scrna_integrated, fid)

# END
