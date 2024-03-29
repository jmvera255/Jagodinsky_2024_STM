#! Rscript

# 03_BayesSpace.R


# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(BayesSpace)
  library(ggplot2)
  library(patchwork)
  library(SpatialExperiment)
})

source("src/BayesSpace.functions.R")

# get cleaned SpatialExperiment object
spe <- qs::qread("BTD3-C1.clean.spe.qs")

# prep spe objects for BayesSpace
spe <- prep_spe(spe)

# viz BayeSpace qTune elbow plot
elbow_plot <- estimate_q(spe)
fid <- "BTD3-C1.BayesSpace.q_elbow_plot.png"
png(fid, width = 1000, height = 800, res = 150)
print(elbow_plot)
invisible(dev.off())

# get BayesSpace clusters w/q = 8
clusters <- get_bayesspace_clusters(spe, 8)

# add clusters to spe colData
colData(spe)$BayesSpace_clusters <- clusters[colnames(spe), 'BayesSpace_clusters']

# save spe, will only have spots with Spot_QC == "keep"
fid <- "BTD3-C1.clusters.spe.qs"
qs::qsave(spe, fid)


# END