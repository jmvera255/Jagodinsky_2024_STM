#! Rscript

# 08_Fig3A_plots.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scLANE)
  
})

source("src/Fig3A.functions.R")
source("src/functions_simulation.R")
source("src/my_plotModels.R")
select <- dplyr::select
filter <- dplyr::filter

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

# cluster plot overlay onto image
cluster_colors <- pals::kelly(n = length(unique(BT_coldata$BayesSpace_clusters)))
names(cluster_colors) <- unique(BT_coldata$BayesSpace_clusters)

p2 <- SpatialDimPlot(spe, group.by = "Clusters", cols = cluster_colors, 
                     crop = FALSE, pt.size.factor = 1) +
  theme(text = element_text(size = 7), 
        legend.key = element_rect(fill = "white"), 
        legend.position = "right", 
        legend.text = element_text(size = 6), 
        legend.key.height = unit(0.275, "cm")) +
  coord_fixed(ratio = image_dims[1]/image_dims[2])

# combine p1 and p2 and save as image
p <- p1 + p2
fid <- "BTD3-C1.Fig3A.plots_1-2.png"
png(fid, height = image_dims[1]*1.1, width = image_dims[2]*2.2, res = 200)
print(p)
invisible(dev.off())

# make violin plots of H2-K1 expression binned by BT dist
df <- select(BT_coldata, BT_dose, barcode)
h2k1_df <- data.frame("H2K1" = spe@assays$SCT@counts['H2-K1',], 
                      barcode = colnames(spe))
df <- left_join(df, h2k1_df, by = "barcode")

# perform wilcoxon rank sum test
high_to_low <- seurat_wilcox_test(spe, x = "high", y = "low")
high_to_low_fdr <- filter(high_to_low, Gene == "H2-K1") %>% 
  pull(p_val_adj)
high_to_low_fdr <- ifelse(high_to_low_fdr[1] < 1E-3, "*", "ns")

high_to_mod <- seurat_wilcox_test(spe, x = "high", y = "low")
high_to_mod_fdr <- filter(high_to_mod, Gene == "H2-K1") %>% 
  pull(p_val_adj)
high_to_mod_fdr <- ifelse(high_to_mod_fdr[1] < 1E-3, "*", "ns")

df$BT_dose <- factor(df$BT_dose, levels = c("high", "mod", "low"))
maxy = ceiling(max(df$H2K1))
p <- ggplot(df, aes(x = BT_dose, y = H2K1)) +
  geom_violin(scale = "width", fill = "grey85") + 
  geom_boxplot(width = 0.15, notch = TRUE) +
  theme_minimal() + ggtitle("BTD3-C1, * FDR < 0.001") +
  ylim(c(0, maxy*1.5)) +
  theme(axis.text = element_text(size = 10)) +
  annotate("segment", y = maxy*1.4, yend = maxy*1.4, 
           x = 1, xend = 3) +
  annotate("text", y = maxy*1.45, x = 2, label = high_to_low_fdr) +
  annotate("segment", y = maxy*1.15, yend = maxy*1.15, 
           x = 1, xend = 2) +
  annotate("text", y = maxy*1.2, x = 1.5, label = high_to_mod_fdr) +
  ylab("H2-K1 expression") + xlab("BT dose region")

fid <- "BTD3-C1.Fig3A.violin.png"
png(fid, height = 750, width = 700, res = 200)
print(p)
invisible(dev.off())

## H2-K1/BT dist plot, Fig 3A plot 4
dist_max <- max(spe$dist_um)
spe[['dist_norm']] <- spe$dist_um/dist_max
order_df <- data.frame(X = spe$dist_norm)
dist_df <- data.frame(X = spe$dist_um)
# there shouldn't be any offset since SCT has normalized for seq depth
cell_offset <- rep(1, ncol(spe))
names(cell_offset) <- colnames(spe)

de_test_glm <- qs::qread("BTD3-C1.scLANE.test_glm.qs")
#de_res_glm <- qs::qread("BTD3-C1.scLANE.res_glm.qs")

p <- my_plotModels(de_test_glm, 
                   gene = "H2-K1", 
                   expr.mat = spe, 
                   plot.gam = FALSE, 
                   plot.glm = FALSE, 
                   plot.null = FALSE,
                   pt = dist_df,
                   size.factor.offset = cell_offset)

y_max <- max(spe@assays$SCT@counts["H2-K1",])

p <- p + xlab(paste0("Distance to BT seed (", "\U03BC", "m)")) +
  ylab("Mhc-1/H2-K1 Expression") + ggtitle("") +
  theme(strip.text = element_blank()) +
  scale_x_continuous(limits = c(0,NA), expand = c(0,0)) +
  ggtitle("BTD3-C1") + 
  scale_y_continuous(limits = c(0, ceiling(y_max*1.3)), 
                     expand = c(0.02, 0), 
                     breaks = seq(0,ceiling(y_max*1.3), by = 10))

# add BT "dose" annotations
p <- p + annotate("segment", y = 1.15*y_max, yend = 1.15*y_max, x = 25, xend = 1985, 
                  color = "black", linewidth = 0.75) +
     annotate("segment", y = 1.15*y_max, yend = 1.15*y_max, x = 2012, xend = 3985, 
              color = "black", linewidth = 0.75) +
     annotate("segment", y = 1.15*y_max, yend = 1.15*y_max, x = 4012, xend = max(dist_df$X), 
              color = "black", linewidth = 0.75) +
     annotate("text", y = 1.22*y_max, x = 600, label = "High Dose", size = 3, hjust = 0) +
     annotate("text", y = 1.22*y_max, x = 2350, label = "Moderate Dose", size = 3, hjust = 0) +
     annotate("text", y = 1.22*y_max, x = 4600, label = "Low Dose", size = 3, hjust = 0)

fid <- "BTD3-C1.Fig3A.scatter.png"
png(fid, height = 1050, width = 1100, res = 250)
print(p)
invisible(dev.off())


# END
