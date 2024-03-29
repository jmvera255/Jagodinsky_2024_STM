#! Rscript

# Fig4D.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated myeloid compartment data
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")

# make Fig5D cell abundance hypergeo heatmap
source("src/cell_abundance_hypergeo_test.R")
ht <- cell_abundance_hypergeo_test(myeloid_integrated, "Myeloid")
fid <- "Fig4D.myeloid.hypergeo_fdr.heatmap.png"
png(fid, width = 8, height = 7.5, units = "cm", res = 300)
draw(ht)
invisible(dev.off())

# END