#! Rscript

# 05_BT_dose_annotations.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(SpatialExperiment)
  library(Seurat)
})

source("src/BT_dose_annotations.functions.R")

# read in SpatialExperiment object which has spatial coordinates stored nicely
# need all spots for calculation thus will use the "raw" spe object created in 02_SpotQC.R
fid <- "BTD3-C1.raw.spe.qs"
spe <- qs::qread(fid)

# pull spatialCoords
spatial_df <- as.data.frame(spatialCoords(spe))
spatial_df$barcode <- rownames(spatial_df)

# get colData
spe_coldata <- as.data.frame(colData(spe))
spe_coldata$barcode <- rownames(spe_coldata)
# read in spatialfactors
fid <- "BTD3-C1/outs/spatial/scalefactors_json.json"
spatial_factors <- jsonlite::read_json(fid)
px_per_um <- spatial_factors$spot_diameter_fullres/65

# load point source annotations
fid <- "../supporting_data/BTD3-C1.pt_source.csv"
pt_source <- read.csv(fid) 
if (nrow(pt_source) == 1) {
  pt_source <- filter(spatial_df, barcode == pt_source$Barcode[1])
} else {
  pt_source <- get_pt_source(spatial_df, spe_coldata, pt_source)
}

# assign labels based on spot distance to point source
test <- define_regions(2000, px_per_um, spatial_df, pt_source)
colnames(test)[grep("test", colnames(test))] <- "high_dose"
test <- define_regions(4000, px_per_um, test, pt_source)
colnames(test)[grep("test", colnames(test))] <- "mod_dose"
test <- define_regions(6000, px_per_um, test, pt_source)
colnames(test)[grep("test", colnames(test))] <- "low_dose"

# do some additional wrangling of dose annotations
in_tissue_spots <- filter(spe_coldata, in_tissue) %>% pull(barcode)
test <- filter(test, barcode %in% in_tissue_spots)
# set default to low dose since some spots may be >6uM from pt source
test$BT_dose <- "low"
test <- mutate(test, 
               BT_dose = ifelse(high_dose, "high", ifelse(mod_dose, "mod", 
                                                          ifelse(low_dose, "low", BT_dose))))

BT_coldata <- left_join(spe_coldata, select(test, barcode, BT_dose), by = "barcode")
BT_coldata$BT_dose[is.na(BT_coldata$BT_dose)] <- "not_tissue"
BT_coldata <- left_join(BT_coldata, spatial_df, by = "barcode")

# add per-spot distance to pt source
BT_coldata$dist_um <- unlist(purrr::map2(BT_coldata$pxl_col_in_fullres, 
                                         BT_coldata$pxl_row_in_fullres, 
                                         function(c, r) calc_distance(c,r, pt_source)))
BT_coldata$dist_um <- BT_coldata$dist_um/px_per_um

# add BayesSpace clusters to BT_coldata
bayesSpace_spe <- qs::qread("BTD3-C1.clusters.spe.qs")
df <- data.frame("BayesSpace_clusters" = bayesSpace_spe$BayesSpace_clusters, 
                 "barcode" = colnames(bayesSpace_spe))
BT_coldata <- left_join(BT_coldata, df, by = "barcode")

# save BT_coldata to file to be used for later analyses
fid <- "BTD3-C1.BT_dose_annotations.qs"
qs::qsave(BT_coldata, fid)

# END