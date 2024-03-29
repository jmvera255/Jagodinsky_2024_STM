#! Rscript

# 04_prep_scRNA_for_TransferData.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(SpatialExperiment)
  library(glmGamPoi)
})

`%notin%` <- Negate(`%in%`)

# get cell cycle marker genes
fid <- "../supporting_data/cell_cycle_genes.qs"
cell_cycle_genes <- qs::qread(fid)

######################################################
### Prep scRNA seq data for transfer of cell types ###
######################################################

# get scrna seurat object with all cell types
fid <- "../scRNA-seq/scrna_qc.scrna_list.processed.qs"
scrna_list <- qs::qread(fid)
sample_names <- names(scrna_list)

scrna_list <- purrr::map(scrna_list, function(x) {
  x <- subset(x, subset = ScDbl_class == "singlet")
  return(x)
})

# remove really low nFeature_RNA cells
scrna_list <- purrr::map(scrna_list, function(x) {
  x <- subset(x, subset = nFeature_RNA >= 350)
  return(x)
})

# remove NA and other cell type prediction cells
scrna_list <- purrr::map(scrna_list, function(scrna) {
  scrna <- subset(scrna, subset = Alt_cell_type != "NA")
  scrna <- subset(scrna, subset = Alt_cell_type != "other")
  return(scrna)
})

# normalize
scrna_list <- purrr::map(scrna_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
  return(x)
})
names(scrna_list) <- sample_names

scrna_list <- purrr::map(scrna_list, function(scrna) {
  N <- length(scrna$orig.ident)
  scrna <- CellCycleScoring(scrna, s.features = unlist(cell_cycle_genes$S), 
                            g2m.features = unlist(cell_cycle_genes$`G2/M`))
  scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score
  scrna <- SCTransform(scrna, vst.flavor = "v2", verbose = FALSE, 
                       vars.to.regress = c("CC.Difference", "percent_mito"), 
                       variable.features.n = 4000, ncells = N)
  
  return(scrna)
})

# save to use as input for Seurat::TransferData
fid <- "scrna_list.TransferData_input.qs"
qs::qsave(scrna_list, fid)



# END