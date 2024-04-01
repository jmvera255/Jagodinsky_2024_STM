
# Fig4F.functions.R

compute_gene_pct <- function(mat, qprob = 0, filter_fun = NULL) {
  gene_quant <- matrixStats::rowQuantiles(mat, probs = qprob)
  
  ncells <- ncol(mat)
  rowSums(mat > gene_quant) / ncells
}

condition_DE <- function(scrna, cluster) {
  # select cells of interest
  scrna <- subset(scrna, subset = seurat_clusters == cluster)
  # subset for variable genes
  scrna <- FindVariableFeatures(scrna, selection.method = "vst", 
                                nfeatures = 7500, assay = "RNA")
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
  results <- tibble(Sample = c("ICI", "BT", "BTICI", "2ICI", "8ICI", "20ICI"))
  results <- mutate(results, Results = purrr::map(Sample, 
                                                  function(x) {DESeq_LRT(counts, scrna$orig.ident, pct_per_sample, test_sample = x)}))
  
  return(results)
}

DESeq_LRT <- function(counts, orig.ident, pct, test_sample) {
  # select for samples of interest
  keep <- c(which(orig.ident == "Sham"), which(orig.ident == test_sample))
  counts <- counts[,keep]
  
  # create sample meta/colData
  cell_meta <- data.frame("Condition" = orig.ident[keep])
  cell_meta <- mutate(cell_meta, Condition = ifelse(Condition == "Sham", "Sham", "Treatment"))
  cell_meta$Condition <- factor(cell_meta$Condition, levels = c("Sham", "Treatment"))
  
  # create dds object and compute sum factors outside DESeq
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = cell_meta,
                                design = ~ Condition)
  qmin_umi = 0.05
  sum_umi <- colSums(SummarizedExperiment::assay(dds))
  min_umi <- quantile(sum_umi, probs = qmin_umi)
  cell_idx <- sum_umi >= min_umi
  dds <- scran::computeSumFactors(dds[,cell_idx])
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1, fitType = "glmGamPoi", 
               minReplicatesForReplace = Inf, useT = TRUE, minmu = 1e-6)
  results <- results(dds, tidy = TRUE, name = "ConditionTreatment")
  results <- select(results, row, log2FoldChange, pvalue, padj)
  colnames(results)[1:2] <- c("Gene", "LFC")
  results <- left_join(results, select(pct, Gene, pct_Sham, 
                                       all_of(paste0("pct_", test_sample))), 
                       by = "Gene")
  idx <- which(colnames(results) %in% c("LFC", "pvalue", "padj"))
  colnames(results)[idx] <- paste(colnames(results)[idx], test_sample, sep = "_")
  
  # for genes with no observed expression in Sham, set LFC to NA
  col_idx <- grep("LFC_", colnames(results))
  row_idx <- which(results$pct_Sham == 0)
  results[row_idx, col_idx] <- NA
  
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

combine_cluster_results <- function(DESeq_LRT) {
  out <- left_join(DESeq_LRT$Results[[1]], 
                   DESeq_LRT$Results[[2]], 
                   by = c("Gene", "pct_Sham"))
  out <- left_join(out, DESeq_LRT$Results[[3]], by = c("Gene", "pct_Sham"))
  out <- left_join(out, DESeq_LRT$Results[[4]], by = c("Gene", "pct_Sham"))
  out <- left_join(out, DESeq_LRT$Results[[5]], by = c("Gene", "pct_Sham"))
  out <- left_join(out, DESeq_LRT$Results[[6]], by = c("Gene", "pct_Sham"))
  
  out <- select(out, Gene, pct_Sham, everything())
  
  return(out)
}

count_hits <- function(df) {
  df <- filter(df, pct_Sham > 0)
  out <- data.frame("hits_ICI" = nrow(filter(df, pct_ICI >= 0.1 & padj_ICI < 0.05)), 
                    "hits_BT" = nrow(filter(df, pct_BT >= 0.1 & padj_BT < 0.05)), 
                    "hits_BTICI" = nrow(filter(df, pct_BTICI >= 0.1 & padj_BTICI < 0.05)), 
                    "hits_2ICI" = nrow(filter(df, pct_2ICI >= 0.1 & padj_2ICI < 0.05)), 
                    "hits_8ICI" = nrow(filter(df, pct_8ICI >= 0.1 & padj_8ICI < 0.05)),
                    "hits_20ICI" = nrow(filter(df, pct_20ICI >= 0.1 & padj_20ICI < 0.05)))
  
  return(out)
}