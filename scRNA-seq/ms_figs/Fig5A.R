#! Rscript

# Fig5A.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated lymphoid compartment data
lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

source("src/make_umap.R")
p <- make_umap(lymphoid_integrated)
fid <- "Fig5A.lymphoid.cell_type.umap.png"
png(fid, height = 900, width = 2000, res = 200)
print(p)
invisible(dev.off())

# END