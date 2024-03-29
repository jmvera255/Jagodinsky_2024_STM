#! Rscript

# Fig4A.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated myeloid compartment data
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")

# make Fig4A UMAP plots
source("src/make_umap.R")
p <- make_umap(myeloid_integrated)
fid <- "Fig4A.myeloid.cell_type.umap.png"
png(fid, height = 900, width = 2000, res = 200)
print(p)
invisible(dev.off())


# END