#! Rscript

# Fig5B.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated lymphoid compartment data
lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

# make Fig5B heatmap
source("src/cluster_Sham_abundance_ht.R")
ht <- cluster_Sham_abundance_ht(lymphoid_integrated, "Lymphoid")
fid <- "Fig5B.cluster_Sham_abundance.lymphoid.heatmap.png"
png(fid, width = 12, height = 10.5 + 4, units = "cm", res = 200)
draw(ht)
invisible(dev.off())

# END