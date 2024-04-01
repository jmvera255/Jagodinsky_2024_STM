#! Rscript

# Fig4E.R

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(glmGamPoi)
  library(DESeq2)
  library(Seurat)
  library(ComplexHeatmap)
  library(ggplot2)
})

######################
##### FUNCTIONS ######
######################

compute_gene_pct <- function(deseq, qprob = 0, filter_fun = NULL) {
  if (!is.null(filter_fun)) {
    deseq <- deseq[, filter_fun(deseq)]
  }
  
  gene_mat <- SummarizedExperiment::assay(deseq)
  gene_quant <- matrixStats::rowQuantiles(gene_mat, probs = qprob)
  
  ncells <- ncol(gene_mat)
  rowSums(gene_mat > gene_quant) / ncells
}

compute_gene_mean <- function(scrna, gene) {
  
  data <- GetAssayData(scrna, assay = "SCT", slot = "data")
  
  gene_quant <- matrixStats::rowQuantiles(gene_mat, probs = qprob)
  
  ncells <- ncol(gene_mat)
  rowSums(gene_mat > gene_quant) / ncells
}

######################
#####    MAIN   ######
######################

# load clustered myeloid compartment object
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")

# select cells of interest
myeloid_integrated <- subset(myeloid_integrated, subset = seurat_clusters == "8")
myeloid_integrated <- subset(myeloid_integrated, subset = orig.ident == "ICI" | orig.ident == "BT" | orig.ident == "BTICI" | orig.ident == "Sham")
invisible(gc())

# get raw RNA assay count matrix for DE analysis
counts <- Seurat::GetAssayData(myeloid_integrated, assay = "RNA", slot = "counts")
counts <- as.matrix(counts)

# subset for variable genes
myeloid_integrated <- FindVariableFeatures(myeloid_integrated, selection.method = "vst", 
                              nfeatures = 4000, assay = "RNA")
counts <- counts[VariableFeatures(myeloid_integrated, assay = "RNA"),]

# create sample meta/colData
cell_meta <- data.frame("Sample" = myeloid_integrated$orig.ident)
cell_meta <- mutate(cell_meta, ICI = ifelse(grepl("ICI", Sample), "1", "0"), 
                    BT = ifelse(grepl("BT", Sample), "1", "0"))
# create dds object and compute sum factors outside DESeq
dds <- DESeqDataSetFromMatrix(countData = counts, colData = cell_meta,
                              design = ~ BT*ICI)
qmin_umi = 0.05
sum_umi <- colSums(SummarizedExperiment::assay(dds))
min_umi <- quantile(sum_umi, probs = qmin_umi)
cell_idx <- sum_umi >= min_umi
dds <- scran::computeSumFactors(dds)
dds <- DESeq(dds, test = "LRT", reduced = ~ BT + ICI, fitType = "glmGamPoi", 
             minReplicatesForReplace = Inf, useT = TRUE, minmu = 1e-6)

pct_Sham <- compute_gene_pct(dds[, dds$Sample == "Sham"], qprob = 0)
pct_ICI <- compute_gene_pct(dds[, dds$Sample == "ICI"], qprob = 0)
pct_BT <- compute_gene_pct(dds[, dds$Sample == "BT"], qprob = 0)
pct_BTICI <- compute_gene_pct(dds[, dds$Sample == "BTICI"], qprob = 0)

results <- tibble::as_tibble(results(dds, tidy = TRUE, name = "BT1.ICI1")) %>%
  dplyr::mutate(pct_Sham, pct_ICI, pct_BT, pct_BTICI)
arrange(results, padj) %>% head()

# Define a set of "hits" 
hits <- filter(results, padj < 0.05 & pct_BTICI > 0.1)

# reload clustered myeloid compartment object to get neutrophil cells from all conditions
myeloid_integrated <- readRDS("myeloid_compartment.integrated.SeuratObj.rds")
# select cells of interest
myeloid_integrated <- subset(myeloid_integrated, subset = seurat_clusters == "8")
myeloid_integrated$orig.ident <- factor(myeloid_integrated$orig.ident, levels = c("Sham", "ICI", "BT", "BTICI", 
                                                        "2ICI", "8ICI", "20ICI"))

# get SCTransformed counts data mat for heatmap visualization
sct_mat <- GetAssayData(myeloid_integrated, assay = "SCT", slot = "data")
sct_mat <- sct_mat[hits$row,]
mean_mat <- as.matrix(cbind(rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "Sham")]), 
                            rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "ICI")]), 
                            rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "BT")]), 
                            rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "BTICI")]), 
                            rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "2ICI")]), 
                            rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "8ICI")]), 
                            rowMeans(sct_mat[,which(myeloid_integrated$orig.ident == "20ICI")])))
rownames(mean_mat) <- hits$row
colnames(mean_mat) <- c("Sham", "ICI", "BT", "BT+ICI", "2 Gy + ICI", 
                        "8 Gy + ICI", "20 Gy + ICI")
gene_means <- rowMeans(mean_mat)
gene_sds <- apply(mean_mat, 1, sd)
gene_zmat <- (mean_mat - gene_means)/gene_sds
rownames(gene_zmat) <- hits$row
colnames(gene_zmat) <- c("Sham", "ICI", "BT", "BT+ICI", "2 Gy + ICI", 
                         "8 Gy + ICI", "20 Gy + ICI")

x <- ceiling(max(abs(c(gene_zmat))))
mat_scale <- seq(-x, x, by = 2)
col_fun <- circlize::colorRamp2(mat_scale, 
                                pals::coolwarm(n = length(mat_scale)),
                                space = "LAB")
ht <- Heatmap(gene_zmat, cluster_rows = TRUE, 
              cluster_columns = TRUE,
              col = col_fun,
              column_title_gp = gpar(fontsize = 10), 
              column_names_gp = gpar(fontsize = 7), 
              row_names_gp = gpar(fontsize = 7), 
              name = "gene\nz-score", 
              height = unit(5.5, "cm"), width = unit(6.5, "cm"),
              column_dend_height = unit(0.5, "cm"), 
              row_dend_width = unit(0.5, "cm"),
              column_title = "BT+ICI Synergistic Gene Hits", 
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7), 
                                          labels_gp = gpar(fontsize = 6), 
                                          grid_height = unit(1, "cm"), 
                                          grid_width = unit(0.25, "cm")))
fid <- "Fig4E.myeloid_neutrophils.BTICI_synergistic_hits.zscore.heatmap.png"
png(fid, width = 9.65, height = 8.3, units = "cm", res = 200)
draw(ht)
invisible(dev.off())


# END