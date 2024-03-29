
# Fig3A.functions.R
seurat_wilcox_test <- function(spe, x, y) {
  spe <- subset(spe, subset = BT_dose == x | BT_dose == y)
  results <- FindMarkers(spe, ident.1 = x, ident.2 = y, 
                         group.by = "BT_dose", assay = "SCT", 
                         slot = "data", test.use = "wilcox", 
                         verbose = FALSE)
  results$Gene <- rownames(results)
  
  return(results)
}
