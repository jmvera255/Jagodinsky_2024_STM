#! Rscript

# Fig4G.R

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
})

# get jmvera modified JASMINE code that is "sourceable" as a function
if (!dir.exists("JASMINE/")){
  system("git clone -b make_source https://github.com/jmvera255/JASMINE.git")
}
source("JASMINE/JASMINE_V1_11October2021.R")

myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")
myeloid_integrated$orig.ident <- factor(myeloid_integrated$orig.ident, levels = c("Sham", "ICI", "BT",
                                                        "BTICI", "2ICI", "8ICI", "20ICI"))
# subset to monocyte clusters of interest
myeloid_integrated <- subset(myeloid_integrated, subset = seurat_clusters %in% c("4", "11", "13"))

gene_counts <- as.matrix(myeloid_integrated@assays$RNA@counts)
cell_meta <- select(myeloid_integrated@meta.data, orig.ident, seurat_clusters)
colnames(cell_meta) <- c("orig.ident", "cluster")
cell_meta$barcode <- colnames(myeloid_integrated)

classical_genes <- scan("../../supporting_data/Orecchioni_2019/classically_activated.txt", 
                        sep = "\n", what = "character")

# Classical activation scoring
scores <- JASMINE(data = gene_counts, 
                  method = 'oddsratio',
                  genes = classical_genes)

scores <- left_join(scores, cell_meta, by = c("SampleID" = "barcode"))
scores$SampleID <- factor(scores$SampleID, 
                          levels = colnames(myeloid_integrated))
scores <- arrange(scores, SampleID)
myeloid_integrated[['classical_activation_score']] <- scores$JAS_Scores

p <- VlnPlot(myeloid_integrated, features = "classical_activation_score", 
             group.by = "orig.ident", split.by = "seurat_clusters", 
             pt.size = 0.1)
p[[1]]$labels$x <- ""
p[[1]]$labels$y <- "Classical Activation Score"
p[[1]]$labels$title <- "Myeloid Clusters 4, 11, 13 (Monocytes)"

fid <- "Fig4G.monocyte.classical_activation.violin.png"
png(fid, width = 1100, height = 600, res = 150)
print(p)
invisible(dev.off())

# END