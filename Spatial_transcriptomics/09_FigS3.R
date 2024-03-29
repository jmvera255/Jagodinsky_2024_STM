#! Rscript

# 09_FigS3.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  
})

# get Seurat spatial data object
fid <- "BTD3-C1.seurat_spatial.qs"
spe <- qs::qread(fid)

# get BT_dose annotation colData
fid <- "BTD3-C1.BT_dose_annotations.qs"
BT_coldata <- qs::qread(fid)

# filter to spots in seurat obj
BT_coldata <- filter(BT_coldata, barcode %in% colnames(spe))
rownames(BT_coldata) <- BT_coldata$barcode

# add BT_coldata to spe
spe[["Clusters"]] <- BT_coldata[colnames(spe), "BayesSpace_clusters"]
spe[['dist_um']] <- BT_coldata[colnames(spe), "dist_um"]
spe[['BT_dose']] <- BT_coldata[colnames(spe), "BT_dose"]
spe$BT_dose <- factor(spe$BT_dose, levels = c("high", "mod", "low"))
spe[['BT_dist_mm']] <- spe$dist_um/1000

# pull tissue section image from spatial seurat obj to add to the plot
image_dims <- dim(spe@images$slice1@image)
slice <- SeuratObject::GetImage(spe, mode = "raster")

# get H&E only image
p1 <- ggplot() + ggpubr::background_image(slice) + 
  suppressWarnings(coord_fixed(ratio = image_dims[1]/image_dims[2]))

# BT Distance (mm) 
p2 <- Seurat::SpatialFeaturePlot(spe, pt.size.factor = 1.05,
                                 features = c("BT_dist_mm"), 
                                 crop = FALSE) + 
        scale_fill_gradientn(limits = c(0, 7), 
                             breaks = seq(0,7, by = 1), 
                             colours = pals::cubicl(n = 8), 
                             name = "BT\ndistance\n(mm)") +
        theme(text = element_text(size = 6), 
              legend.direction = "vertical", 
              legend.position = "right", 
              legend.key.width = unit(10, "pt"), 
              legend.key.height = unit(0.33, "cm")) +
        suppressWarnings(coord_fixed(ratio = image_dims[1]/image_dims[2]))

# HT-K2 SpatialFeature plot
p3 <- SpatialFeaturePlot(spe, features = "H2-K1", 
                         crop = FALSE, pt.size.factor = 1) +
        theme(text = element_text(size = 7), 
              legend.direction = "vertical", 
              legend.position = "right", 
              legend.key.width = unit(10, "pt"), 
              legend.key.height = unit(0.3, "cm")) +
        coord_fixed(ratio = image_dims[1]/image_dims[2]) #+
#scale_fill_gradientn(limits = c(0, 4.5), 
#                     breaks = seq(0,4.5, by = 1), 
#                     colours = rev(pals::brewer.spectral(n = 5)), 
#                     name = "H2-K1")


# BT_dose SpatialDim plot
dose_colors <- c("high" = "salmon", "mod" = "#00BA38", "low" = "#619CFF")
p4 <- SpatialDimPlot(spe, group.by = "BT_dose", 
                     cols = dose_colors, 
                     crop = FALSE, pt.size.factor = 1) +
  theme(text = element_text(size = 7), 
        legend.key = element_rect(fill = "white"), 
        legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.key.height = unit(0.275, "cm")) +
  coord_fixed(ratio = image_dims[1]/image_dims[2]) +
  guides(guide_legend(override.aes = list(size = 1.5)))

# combine plots into single patchwork object
p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
fid <- "BTD3-C1.FigS3.png"
png(fid, height = image_dims[1]*2.1, width = image_dims[2]*2.2, res = 200)
print(p)
invisible(dev.off())

# END
