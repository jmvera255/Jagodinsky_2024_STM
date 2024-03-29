#! Rscript

# Fig5C.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(qs)
  library(purrr)
  library(ggplot2)
})

# read in integrated lymphoid compartment data
lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

# make Fig5C cell type dot plot
source("src/cell_type_dot_plot.R")
p <- cell_type_dot_plot(lymphoid_integrated)
fid <- "Fig5C.lymphoid.cell_type.dotplot.png"
png(fid, height = 650, width = 1200)
print(p)
invisible(dev.off())

# END