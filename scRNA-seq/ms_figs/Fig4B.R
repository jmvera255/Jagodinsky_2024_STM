#! Rscript

# Fig4B.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated myeloid compartment data
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")

# make Fig4B heatmap
source("src/cluster_Sham_abundance_ht.R")
ht <- cluster_Sham_abundance_ht(myeloid_integrated, "Myeloid")
fid <- "Fig4B.cluster_Sham_abundance.myeloid.heatmap.png"
png(fid, width = 12, height = 8.5 + 4, units = "cm", res = 200)
draw(ht)
invisible(dev.off())

# END