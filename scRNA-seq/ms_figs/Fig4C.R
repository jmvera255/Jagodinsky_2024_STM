#! Rscript

# Fig4C.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated myeloid compartment data
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")

# make Fig4C cell type dot plot
source("src/cell_type_dot_plot.R")
p <- cell_type_dot_plot(myeloid_integrated)
fid <- "Fig4C.myeloid.cell_type.dotplot.png"
png(fid, height = 650, width = 1200)
print(p)
invisible(dev.off())

# END