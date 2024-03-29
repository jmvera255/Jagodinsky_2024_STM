
suppressPackageStartupMessages({
  library(ggplot2)
})
# make_umap.R
make_umap <- function(scrna) {
  # factor alt cell type label with levels sorted by desc character length
  df <- data.frame(label = unique(scrna$Alt_cell_type))
  df$n <- unlist(purrr::map(df$label, nchar))
  df <- arrange(df, desc(n))
  
  # plot colors
  type_cols <- c(pals::alphabet2(), pals::glasbey())
  names(type_cols) <- NULL
  type_cols <- type_cols[1:length(unique(scrna$Alt_cell_type))]
  names(type_cols) <- unique(sort(scrna$Alt_cell_type))
  type_cols <- type_cols[df$label]
  
  scrna$Alt_cell_type <- factor(scrna$Alt_cell_type, levels = df$label)
  
  idx <- which(scrna$orig.ident %in% c("2ICI", "8ICI", "20ICI", "BTICI"))
  scrna <- scrna[,idx]
  scrna$orig.ident <- stringr::str_replace(scrna$orig.ident, "BTICI", "BT+ICI")
  scrna$orig.ident <- stringr::str_replace(scrna$orig.ident, "2ICI", "2 Gy + ICI")
  scrna$orig.ident <- stringr::str_replace(scrna$orig.ident, "8ICI", "8 Gy + ICI")
  scrna$orig.ident <- stringr::str_replace(scrna$orig.ident, "20ICI", "20 Gy + ICI")
  scrna$orig.ident <- factor(scrna$orig.ident, 
                             levels = c("BT+ICI", "2 Gy + ICI", "8 Gy + ICI", "20 Gy + ICI"))
  
  # make fig4 panel A umap plot
  p <- DimPlot(scrna, reduction = "umap", split.by = "orig.ident",
               group.by = "Alt_cell_type", cols = alpha(type_cols, 0.7),
               ncol = 4, 
               pt.size = 0.35)
  p$labels$title <- NULL
  p <- p + theme(axis.text = element_text(size = 8.5),
                 axis.title = element_text(size = 8.5),
                 title = element_text(size = 8), 
                 legend.text = element_text(size = 8), 
                 legend.position = "bottom") +
    guides(col = guide_legend(ncol = 5, override.aes = list(size = 1), 
                              keyheight = 0.65))
  return(p)
}