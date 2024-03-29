
# SpotQC.functions.R

calc_background <- function(spe) {
  df <- as.data.frame(table(spe$in_tissue)) %>% 
    tidyr::pivot_wider(names_from = Var1, values_from = Freq)
  colnames(df) <- c("non_tissue", "tissue")
  df$sum <- length(spe$in_tissue)
  df$pct_background <- df$non_tissue[1]/df$sum[1]
  
  return(df)
}

annotate_coldata <- function(spe, sample) {
  i <- grep(sample, sample_filters$sample_id)
  
  colData(spe)$nFeatures_group <- 
    ifelse(colData(spe)$nFeatures < sample_filters$nFeatures_threshold[i], "drop", "keep")
  colData(spe)$nCounts_group <- 
    ifelse(colData(spe)$decont_nCount > sample_filters$nCounts_threshold[i], "drop", "keep")
  colData(spe)$manual_QC <- "keep"
  
  return(spe)
}
