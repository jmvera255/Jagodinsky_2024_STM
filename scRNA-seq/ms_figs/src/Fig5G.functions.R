
# Fig5G.functions.R

compute_gene_pct <- function(mat, qprob = 0, filter_fun = NULL) {
  gene_quant <- matrixStats::rowQuantiles(mat, probs = qprob)
  
  ncells <- ncol(mat)
  rowSums(mat > gene_quant) / ncells
}

synergistic_DE <- function(scrna, cluster) {
  # select cells of interest
  scrna <- subset(scrna, subset = seurat_clusters == cluster)
  # subset for variable genes
  scrna <- FindVariableFeatures(scrna, selection.method = "vst", 
                                nfeatures = 5000, assay = "RNA")
  # make counts matrix with all samples
  counts <- GetAssayData(scrna, assay = "RNA", slot = "counts")
  counts <- as.matrix(counts)
  counts <- counts[VariableFeatures(scrna, assay = "RNA"),]
  
  # make df of per gene sample percentages
  pct_per_sample <- data.frame("Gene" = rownames(counts), 
                               "pct_Sham" = compute_gene_pct(counts[,which(scrna$orig.ident == "Sham")]),
                               "pct_ICI" = compute_gene_pct(counts[,which(scrna$orig.ident == "ICI")]),
                               "pct_BT" = compute_gene_pct(counts[,which(scrna$orig.ident == "BT")]),
                               "pct_BTICI" = compute_gene_pct(counts[,which(scrna$orig.ident == "BTICI")]),
                               "pct_2ICI" = compute_gene_pct(counts[,which(scrna$orig.ident == "2ICI")]),
                               "pct_8ICI" = compute_gene_pct(counts[,which(scrna$orig.ident == "8ICI")]),
                               "pct_20ICI" = compute_gene_pct(counts[,which(scrna$orig.ident == "20ICI")]))
  results <- tibble(Sample = c("BTICI", "2ICI", "8ICI", "20ICI"))
  results <- mutate(results, Results = purrr::map(Sample, 
                                                  function(x) {DESeq_LRT(counts, scrna$orig.ident, pct_per_sample, test_sample = x)}))
  
  return(results)
}

DESeq_LRT <- function(counts, orig.ident, pct, test_sample) {
  # select for samples of interest
  keep <- c(which(orig.ident == "Sham"), which(orig.ident == "ICI"), 
            which(orig.ident == "BT"), which(orig.ident == test_sample))
  counts <- counts[,keep]
  
  # create sample meta/colData
  cell_meta <- data.frame("Sample" = orig.ident[keep])
  cell_meta <- mutate(cell_meta, ICI = ifelse(grepl("ICI", Sample), "1", "0"), 
                      RAD = ifelse(grepl("BT|8|2", Sample), "1", "0"))
  cell_meta$RAD <- factor(cell_meta$RAD)
  cell_meta$ICI <- factor(cell_meta$ICI)
  
  # create dds object and compute sum factors outside DESeq
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = cell_meta,
                                design = ~ RAD*ICI)
  qmin_umi = 0.05
  sum_umi <- colSums(SummarizedExperiment::assay(dds))
  min_umi <- quantile(sum_umi, probs = qmin_umi)
  cell_idx <- sum_umi >= min_umi
  dds <- scran::computeSumFactors(dds)
  dds <- DESeq(dds, test = "LRT", reduced = ~ RAD + ICI, fitType = "glmGamPoi", 
               minReplicatesForReplace = Inf, useT = TRUE, minmu = 1e-6)
  results <- results(dds, tidy = TRUE, name = "RAD1.ICI1")
  results <- left_join(results, select(pct, Gene, pct_Sham, pct_ICI, pct_BT, matches(test_sample)), 
                       by = c("row" = "Gene"))
  results <- select(results, -lfcSE)
  
  return(results)
}

process_results <- function(cluster, compartment, df, scrna) {
  hits <- purrr::map(df$Results, get_hits)
  plot_results(cluster, compartment, hits, scrna)
}

get_hits <- function(df) {
  col_str <- colnames(df)[ncol(df)]
  colnames(df)[ncol(df)] <- "test"
  out <- filter(df, padj < 0.05 & test > 0.1)# %>% pull(row)
  colnames(df)[ncol(df)] <- col_str
  
  return(out)
}

plot_results <- function(cluster, compartment, hits, scrna) {
  genes <- unlist(purrr::map(hits, function(df) {return(df$row)})) %>% unique()
  scrna <- subset(scrna, subset = seurat_clusters == cluster)
  # get hits that are also in SCT assay
  genes <- table(c(genes, rownames(GetAssayData(scrna, assay = "SCT", slot = "counts")))) %>%
    as.data.frame() %>% filter(Freq == 2) %>% pull(Var1) %>% as.character()
  # filter scrna obj for genes of interest
  scrna <- scrna[genes,]
  
  mat <- gene_zmat[genes,]
  h =  nrow(mat)/5
  
  bot_ha <- columnAnnotation(Condition = scrna$orig.ident, 
                             col = list(Condition = c("Sham" = sample_colors[2], 
                                                      "ICI" = sample_colors[3], 
                                                      "BT" = sample_colors[4], 
                                                      "BTICI" = sample_colors[5], 
                                                      "2ICI" = sample_colors[6], 
                                                      "8ICI" = sample_colors[7], 
                                                      "20ICI" = sample_colors[8])))
  ht <- Heatmap(mat, cluster_rows = TRUE, 
                cluster_column_slices = TRUE,
                column_split = scrna$orig.ident,
                name = "gene z-score",
                column_title = paste(compartment, "cluster", cluster, 
                                     "Interaction Term Genes"), 
                row_title = "Genes", 
                row_title_gp = gpar(fontsize = 14), 
                height = unit(h, units = "cm"), 
                width = unit(11, units = "cm"), 
                row_names_gp = gpar(fontsize = 8), 
                bottom_annotation = bot_ha, 
                show_column_names = FALSE)
  fid <- paste0("pngs/", fid.prefix, ".", compartment, ".cluster_", 
                cluster, ".z-score_heatmap.png")
  png(fid, height = h + 5, width = 15, units = "cm", res = 200)
  draw(ht)
  invisible(dev.off())
}

annotate_EBRT_padj <- function(df) {
  out <- df$Results[[which(df$Sample == "BTICI")]]
  colnames(out)[grep("padj", colnames(out))] <- "padj_BTICI"
  colnames(out)[grep("row", colnames(out))] <- "Gene"
  
  # 2ICI
  EBRT_df <- select(df$Results[[which(df$Sample == "2ICI")]], padj, pct_2ICI)
  out$`padj_2ICI` <- EBRT_df$padj
  out$pct_2ICI <- EBRT_df$pct_2ICI
  
  # 8ICI
  EBRT_df <- select(df$Results[[which(df$Sample == "8ICI")]], padj, pct_8ICI)
  out$`padj_8ICI` <- EBRT_df$padj
  out$pct_8ICI <- EBRT_df$pct_8ICI
  
  # 20ICI
  EBRT_df <- select(df$Results[[which(df$Sample == "20ICI")]], padj, pct_20ICI)
  out$`padj_20ICI` <- EBRT_df$padj
  out$pct_20ICI <- EBRT_df$pct_20ICI
  
  return(out)
}

count_hits <- function(df) {
  df <- filter(df, padj_BTICI < 0.05 & pct_BTICI > 0.1)
  out <- data.frame("hits_BTICI" = nrow(df), 
                    "hits_2ICI" = length(which(df$padj_2ICI < 0.05)), 
                    "hits_8ICI" = length(which(df$padj_8ICI < 0.05)), 
                    "hits_20ICI" = length(which(df$padj_20ICI < 0.05)))
  df <- filter(df, padj_2ICI >= 0.05 & padj_8ICI >= 0.05 & padj_20ICI >= 0.05)
  out$unique_BTICI <- nrow(df)
  
  return(out)
}

plot_interaction_hits <- function(cluster, results, scrna, hits = NULL) {
  DefaultAssay(scrna) <- "SCT"
  scrna <- subset(scrna, subset = seurat_clusters == cluster)
  title <- paste0("Lymphoid Cluster ", cluster, 
                  "\nInteraction Term Genes")
  
  if (is.null(hits)) {
    hits <- get_hits(results$Results[[1]])
    hits <- hits$row
    title <- paste0("Lymphoid Cluster ", cluster, 
                    "\nBTICI Interaction Term Genes")
  }
  both <- table(c(rownames(scrna), hits)) %>% 
    as.data.frame() %>% filter(Freq == 2) %>% pull(Var1) %>% as.character()
  scrna <- scrna[both,]
  
  mat <- get_zscore_mat(scrna)
  h <- nrow(mat)/4.5
  
  x <- ceiling(max(abs(c(mat))))
  mat_scale <- seq(-x, x, by = 2)
  col_fun <- circlize::colorRamp2(mat_scale, 
                                  pals::coolwarm(n = length(mat_scale)),
                                  space = "LAB")
  ht <- Heatmap(mat, cluster_rows = TRUE, 
                cluster_columns = TRUE, 
                name = "gene z-score", 
                col = col_fun,
                row_names_side = "left", 
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 10),
                column_title_gp = gpar(fontsize = 12),
                column_title = title, 
                height = unit(h, "cm"), width = unit(5, "cm"))

  return(ht)
}

get_zscore_mat <- function(scrna) {
  # Prep SCtrans count data for heatmap plotting
  sct_mat <- GetAssayData(scrna, assay = "SCT", slot = "data")
  sct_mat <- as.matrix(sct_mat)
  mean_mat <- as.matrix(cbind(rowMeans(sct_mat[,which(scrna$orig.ident == "Sham")], na.rm = TRUE), 
                              rowMeans(sct_mat[,which(scrna$orig.ident == "ICI")], na.rm = TRUE), 
                              rowMeans(sct_mat[,which(scrna$orig.ident == "BT")], na.rm = TRUE), 
                              rowMeans(sct_mat[,which(scrna$orig.ident == "BTICI")], na.rm = TRUE), 
                              rowMeans(sct_mat[,which(scrna$orig.ident == "2ICI")], na.rm = TRUE), 
                              rowMeans(sct_mat[,which(scrna$orig.ident == "8ICI")], na.rm = TRUE), 
                              rowMeans(sct_mat[,which(scrna$orig.ident == "20ICI")], na.rm = TRUE)))
  rownames(mean_mat) <- rownames(sct_mat)
  colnames(mean_mat) <- c("Sham", "ICI", "BT", "BTICI", "2ICI", "8ICI", "20ICI")
  gene_means <- rowMeans(mean_mat, na.rm = TRUE)
  gene_sds <- apply(mean_mat, 1, sd, na.rm = TRUE)
  gene_zmat <- (mean_mat - gene_means)/gene_sds
  rownames(gene_zmat) <- rownames(sct_mat)
  colnames(gene_zmat) <- c("Sham", "ICI", "BT", "BTICI", "2ICI", "8ICI", "20ICI")
  
  return(gene_zmat)
}