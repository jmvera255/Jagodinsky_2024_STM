#! Rscript

# Fig5D.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated lymphoid compartment data
lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

# make Fig5D cell abundance hypergeo heatmap
source("src/cell_abundance_hypergeo_test.R")
ht <- cell_abundance_hypergeo_test(lymphoid_integrated, "Lymphoid")
fid <- "Fig5D.lymphoid.hypergeo_fdr.heatmap.png"
png(fid, width = 8, height = 7.75, units = "cm", res = 300)
draw(ht)
invisible(dev.off())


# END