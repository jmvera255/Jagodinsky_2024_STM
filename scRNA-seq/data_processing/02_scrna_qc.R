#! /usr/bin/env bash
# 02_scrna_qc.R

# add functions needed for this work
source('src/scrna_qc.functions.R')

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(scDblFinder)
  library(SingleR)
})

# define vector of each sample name/label
samples <- c("Sham_S1", "BT_S2","BTICI_S3", "ICI_S4", "2ICI_S5", "8ICI_S6", "20ICI_S7")

# load filtered count matrix data from cell ranger to build seurat objects
data10x <- list()
for (s in samples) {
  dir <- paste0(s, "/outs/filtered_feature_bc_matrix/")
  data10x[[s]] <- Read10X(dir)
}

# remove generic sample id from sample labels
samples <- stringr::str_split(samples, "_", simplify = TRUE)[,1]
names(data10x) <- samples
scrna_list <- purrr::map(samples, function(x) CreateSeuratObject(counts = data10x[[x]], 
                                                                 project = x,
                                                                 min.cells = 10, 
                                                                 min.features = 100))
names(scrna_list) <- samples
rm(data10x)

# calc percent mito and remove >10%
scrna_list <- purrr::map(scrna_list, function(x) {
  x[["percent_mito"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x <- subset(x, subset = percent_mito < 10)
  return(x)
})

# convert Seurat object to sce to run scDblFinder
sce_list <- purrr::map(scrna_list, as.SingleCellExperiment)
sce_list <- purrr::map(sce_list, function(x) {
  # set.seed() # use set.seed for consistent doublet prediction run-to-run
  x <- scDblFinder(x)
  return(x)
})

# mark singlet/doublet predictions in seurat object list
scrna_list <- purrr::map2(scrna_list, sce_list, function(x, y) {
  x[['ScDbl_class']] <- y$scDblFinder.class
  return(x)
})

# perform cell type predictions with singleR using mouse ImmGen reference from celldex
cell_onto <- ontoProc::getCellOnto(useNew = TRUE)
cell_onto_df <- data.frame("CL" = names(cell_onto$name), 
                           "Cell_type" = cell_onto$name)

mouse_imm_ref <- celldex::ImmGenData(cell.ont = "all")
mouse_imm_df <- data.frame(CL = mouse_imm_ref$label.ont, 
                           Label = mouse_imm_ref$label.fine,
                           Name = mouse_imm_ref$label.main)
mouse_imm_df <- unique(mouse_imm_df)
mouse_imm_df <- left_join(mouse_imm_df, cell_onto_df, by = "CL")

singleR_annotations <- purrr::map(scrna_list, function(x) {
                        run_singleR(x, mouse_imm_ref, "mouse_imm_ref")})

# do initial addition of singleR predictions
scrna_list <- purrr::map2(scrna_list, singleR_annotations, process_cell_annotations)

# rework the singleR annotations
singleR_annotations <- purrr::map2(seq_along(singleR_annotations), singleR_annotations, 
                                   function(i, y) {
                                     df <- data.frame(UMI = rownames(y), 
                                                      Sample = i, 
                                                      Label = y$labels, 
                                                      Pruned_label = y$pruned.labels)
                                     return(df)
                                   })
singleR_annotations <- do.call(rbind.data.frame, singleR_annotations)
singleR_annotations <- left_join(singleR_annotations, mouse_imm_df, by = "Label")

# create sample-specific UMIs in singleR_annotations df
singleR_annotations$UMI <- paste(singleR_annotations$UMI, singleR_annotations$Sample, sep = "_")
rownames(singleR_annotations) <- singleR_annotations$UMI

# to improve the cell type visualization I want to provide an alternative label for 
# the cell types that occur less than 10 times in the data set
cell_type_tbl <- table(singleR_annotations$Cell_type) %>% as.data.frame()
cells <- filter(cell_type_tbl, Freq < 10) %>% pull(Var1)
cells <- paste(cells, collapse = "|")
singleR_annotations$Alt_cell_type <- stringr::str_replace_all(singleR_annotations$Cell_type, 
                                                              cells, "other")

# shorten cell type labels
singleR_annotations$Cell_type <- stringr::str_replace_all(singleR_annotations$Cell_type, 
                                                          "-positive", "(+)")
singleR_annotations$Cell_type <- stringr::str_replace_all(singleR_annotations$Cell_type, 
                                                          "-negative", "(-)")
singleR_annotations$Alt_cell_type <- stringr::str_replace_all(singleR_annotations$Alt_cell_type, 
                    "-positive", "(+)")
singleR_annotations$Alt_cell_type <- stringr::str_replace_all(singleR_annotations$Alt_cell_type, 
                                                              "-negative", "(-)")

# remove ScDblFinder predicted doublets
scrna_list <- purrr::map(scrna_list, function(x) {
  x <- subset(x, subset = ScDbl_class == "singlet")
  return(x)
})

# remove cells with nFeature_RNA < 350
scrna_list <- purrr::map(scrna_list, function(x) {
  x <- subset(x, subset = nFeature_RNA >= 350)
  return(x)
})

scrna_list <- purrr::map(scrna_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
  return(x)
})

# add singleR non-pruned labels and other additional modified cell labels
scrna_list <- purrr::map(seq_along(scrna_list), set_new_labels)
names(scrna_list) <- samples

# then save scrna_list with all cell types to use as input for subsequent work
fid <- "scrna_qc.scrna_list.processed.qs"
qs::qsave(scrna_list, fid) #2022-06-13_integration_optimization_prep 

# pull out predicted lymphoid cells
lymphoid_scrna_list <- purrr::map(scrna_list, function(x) {
  x <- subset(x, subset = basic_cell_name == "B cells" | basic_cell_name == "NK cells" |   basic_cell_name == "NKT" | basic_cell_name == "Tgd" | basic_cell_name == "ILC" | basic_cell_name ==   "T cells")
  return(x)
})
names(lymphoid_scrna_list) <- samples
fid <- "scrna_qc.lymphoid.scrna_list.processed.qs"
qs::qsave(lymphoid_scrna_list, fid) #2022-09-02_lymphoid_compartment_clustering scrna_list

# pull out predicted myeloid cells
myeloid_scrna_list <- purrr::map(scrna_list, function(x) {
  x <- subset(x, subset = basic_cell_name == "DC" | basic_cell_name == "Macrophages"   | basic_cell_name == "Mast cells" | basic_cell_name == "Monocytes" | basic_cell_name == "Neutrophils" |   basic_cell_name == "Basophils" | basic_cell_name == "Microglia")
  return(x)
})
names(myeloid_scrna_list) <- samples
fid <- "scrna_qc.myeloid.scrna_list.processed.qs"
qs::qsave(myeloid_scrna_list, fid) #2022-09-06_myeloid_compartment_clustering scrna_list


# END