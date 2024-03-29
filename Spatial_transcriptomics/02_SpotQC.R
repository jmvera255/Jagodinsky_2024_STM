#! Rscript

# 02_SpotQC.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(SpotClean)
  library(SpatialExperiment)
})

source("src/SpotQC.functions.R")

# read in raw data
data_dir <- paste0("SpaceRanger/BTD3-C1/outs/")
spe <- read10xVisium(samples = data_dir, sample_id = "BTD3-C1", 
                     type = "HDF5", images = "hires", data = "raw")

# remove deprecated ensembl gene entries
idx <- grep("DEPRECATED", rownames(assays(spe)$counts), invert = TRUE)
spe <- spe[idx,]

# remove reliance on H5FD file and allow spe object to be portable
assays(spe)$counts <- as(assays(spe)$counts, "CsparseMatrix")

# collect some per spot metrics
df <- data.frame(in_tissue = spe$in_tissue, 
                 nCount = colSums(assays(spe)$counts))
df$raw_log2_nCount <- log2(df$nCount + 1)

# add to colData slot
colData(spe)$raw_log2_nCount <- df$raw_log2_nCount
colData(spe)$raw_nCount <- df$raw_nCount

# add to metadata slot
metadata(spe)$slide$raw_log2_nCount <- df$raw_log2_nCount
metadata(spe)$slide$raw_nCount <- df$raw_nCount

# Calc background
tissue_spot_metrics <- calc_background(spe)
tissue_spot_metrics$Sample <- "BTD3-C1"
metadata(spe)$pct_nontissue <- tissue_spot_metrics$pct_background[1]

# test for at least 25% non-tissue/background spots
if (tissue_spot_metrics$pct_background[1] < 0.25) {
  
  # save raw spe object
  assays(spe)$counts <- as(assays(spe)$counts, "CsparseMatrix")
  fid <- "BTD3-C1.raw.spe.qs"
  qs::qsave(spe, fid)
  
  stop(paste0("The percent of non-tissue spots is less than 25%\n",
              "and SpotClean is not recommended in this case."))
}

# Run SpotClean 
spe_clean <- spotclean(spe, gene_cutoff = 0)

# save raw Spatial Experiment Object
assays(spe)$counts <- as(assays(spe)$counts, "CsparseMatrix")
fid <- "BTD3-C1.raw.spe.qs"
qs::qsave(spe, fid)

# collect some spe_clean per spot metrics
df <- data.frame(in_tissue = spe_clean$in_tissue, 
                 decont_nCount = colSums(assays(spe_clean)$decont))
df$decont_log2_nCount <- log2(df$decont_nCount + 1)

# add to colData slot
colData(spe_clean)$decont_log2_nCount <- df$decont_log2_nCount
colData(spe_clean)$decont_nCount <- df$decont_nCount

# add to metadata slot
metadata(spe_clean)$slide$decont_log2_nCount <- df$decont_log2_nCount
metadata(spe_clean)$slide$decont_nCount <- df$decont_nCount

## Per spot QC

# load predefined metric thresholds
sample_filters <- read.delim("../supporting_data/ST_spotQC_thresholds.txt")

# calculate nFeatures for each spot
colData(spe_clean)$nFeatures <- colSums(assays(spe_clean)$decont[,] > 0)

# flag low qual spots, creates manual_QC colData
spe_clean <- annotate_coldata(spe_clean, "BTD3-C1")

# flag manually identified spots
for (umi in c("GTGTGAGCCGAGGTGC-1","CGAGATGTTGCCTATA-1", 
              "TATGGAGTTTCTCGTT-1", "CTTGTCGTACGTGTCA-1")) {
  i <- grep(umi, rownames(colData(spe_clean)))
  spe_clean$manual_QC[i] <- "drop"
}

# save cleaned spe object
assays(spe_clean)$counts <- as(assays(spe_clean)$counts, "CsparseMatrix")
fid <- "BTD3-C1.clean.spe.qs"
qs::qsave(spe_clean, fid)

# END
