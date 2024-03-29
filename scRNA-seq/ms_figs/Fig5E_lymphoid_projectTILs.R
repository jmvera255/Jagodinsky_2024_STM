#! Rscript

# Fig5E_lymphoid_projectTILs.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ProjecTILs)
  library(ggplot2)
})

# load mouse TIL Atlas v1 to use as projection reference
download.file(url = "https://doi.org/10.6084/m9.figshare.12478571", 
              destfile = "ref_TILAtlas_mouse_v1.rds")
ref <- readRDS("ref_TILAtlas_mouse_v1.rds")

# read in integrated lymphoid compartment data
lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

# split into list object
scrna_list <- SplitObject(lymphoid_integrated, split.by = "orig.ident")

# make projections
scrna_list_projected <- purrr::map(scrna_list, function(x) {
  make.projection(x, ref = ref)
})

# add cell type predictions
scrna_list_projected <- purrr::map(scrna_list_projected, function(x) {
  cellstate.predict(ref = ref, query = x)
})

# combine all the TIL predictions into a df for use in other analyses
projectTILs_df <- purrr::map(names(scrna_list_projected), function(sample) {
                    return(data.frame(UMI = names(scrna_list_projected[[sample]]$orig.ident),
                                      Sample = sample, 
                                      TIL = scrna_list_projected[[sample]]$functional.cluster))
                  })
projectTILs_df <- do.call(rbind.data.frame, projectTILs_df)

sample_key <- 1:7
names(sample_key) <- c("Sham", "BT", "BTICI", "ICI", "2ICI", "8ICI", "20ICI")
projectTILs_df$I <- recode(projectTILs_df$Sample, !!!sample_key)
projectTILs_df <- mutate(projectTILs_df, UMI = paste(UMI, I, sep = "_"))

# now load integrated lymphoid data
scrna_integrated <- qs::qread("lymphoid_integrated.scrna.qs")

# make projectTIL dot plot Fig5E
scrna_df <- data.frame(UMI = names(scrna_integrated$orig.ident), 
                       Sample = scrna_integrated$orig.ident, 
                       Cluster = scrna_integrated$seurat_clusters)
scrna_df$I <- recode(scrna_df$Sample, !!!sample_key)
scrna_df <- mutate(scrna_df, UMI = paste(UMI, I, sep="_"))
scrna_df <- left_join(projectTILs_df, dplyr::select(scrna_df, -I, -Sample), by = "UMI")
cluster_sums <- as.data.frame(table(scrna_df$Cluster))

scrna_df <- group_by(scrna_df, Cluster, TIL) %>% summarise(N = n(), .groups = "drop")
scrna_df <- left_join(scrna_df, cluster_sums, by = c("Cluster" = "Var1"))

colnames(scrna_df) <- c("cluster", "cell_type", "cell_count", "cluster_total_cells")
scrna_df$pct_type <- scrna_df$cell_count/scrna_df$cluster_total_cells

p <- ggplot(scrna_df, aes(y = cluster, x = cell_type, 
                          size = pct_type*100, color = pct_type*100)) +
  geom_point() +
  scale_color_gradient(high = "red2", low = "grey80") +
  labs(title = "TILPRED abundance by cluster", 
       y = "Lymphoid Cluster", x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   hjust = 1, size = 9), 
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8), 
        title = element_text(size = 12, hjust = 0.5))
p$labels$size <- "% of\ncluster"
p$labels$colour <- "% of\ncluster"

fid <- "Fig5E.lymphoid.TILPRED.dot_plot.png"
png(fid, height = 1150, width = 1100, res = 200)
print(p)
invisible(dev.off())


# END
