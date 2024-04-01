#! Rscript

# Fig5F.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(Seurat)
  library(ggplot2)
})

# get jmvera modified JASMINE code that is "sourceable" as a function
if (!dir.exists("JASMINE/")){
  system("git clone -b make_source https://github.com/jmvera255/JASMINE.git")
}
source("JASMINE/JASMINE_V1_11October2021.R")

# get integrated lymphoid compartment data
lymphoid_integrated <- readRDS("lymphoid_compartment.integrated.SeuratObj.rds")

# get KEGG gene sets
gene_map <- qs::qread("../../supporting_data/Fig4H.gene_map.qs")
gene_map <- select(gene_map, -Symbol)

mmu_paths <- qs::qread("../../supporting_data/mmu_paths.qs")
test_paths <- openxlsx::read.xlsx("../../supporting_data/Fig4H.pathways.xlsx")
test_paths <- filter(test_paths, !is.na(`Pathway.(Orig.results)`))
mmu_paths$all$KEGGPATHID2NAME$to <- stringr::str_remove(mmu_paths$all$KEGGPATHID2NAME$to, 
                                                        " - Mus musculus \\(house mouse\\)")
mmu_paths$all$KEGGPATHID2NAME <- filter(mmu_paths$all$KEGGPATHID2NAME, 
                                        to %in% test_paths$`Pathway.(Orig.results)`)
mmu_paths$all$KEGGPATHID2EXTID <- inner_join(mmu_paths$all$KEGGPATHID2EXTID, 
                                             gene_map, by = c("to" = "gene_id"))
mmu_paths$all$KEGGPATHID2EXTID <- filter(mmu_paths$all$KEGGPATHID2EXTID, 
                                         from %in% mmu_paths$all$KEGGPATHID2NAME$from)
gene_sets <- group_by(mmu_paths$all$KEGGPATHID2EXTID, from) %>% tidyr::nest()
gene_set_dict <- mmu_paths$all$KEGGPATHID2NAME$to
names(gene_set_dict) <- mmu_paths$all$KEGGPATHID2NAME$from

gene_counts <- as.matrix(lymphoid_integrated@assays$RNA@counts)

cell_meta <- select(lymphoid_integrated@meta.data, orig.ident, seurat_clusters)
cell_meta$SampleID <- colnames(lymphoid_integrated)

results <- purrr::map2(gene_sets$from, gene_sets$data, function(id, gs, meta = cell_meta) {
  out <- JASMINE(data = gene_counts, method = 'oddsratio',
                 genes = gs$symbol)
  out$gs_id <- id
  out <- left_join(out, meta, by = "SampleID")
  return(out)
})


results <- do.call(rbind.data.frame, results)

results_sum <- group_by(results, gs_id, orig.ident) %>% 
  summarize(mean_score = mean(JAS_Scores), .groups = "drop")

results_sum <- tidyr::pivot_wider(results_sum, names_from = gs_id, values_from = mean_score)
mat <- as.matrix(t(select(results_sum, -orig.ident)))
colnames(mat) <- results_sum$orig.ident
rownames(mat) <- gene_set_dict[colnames(results_sum)[-1]]

# define heatmap color scale
mat_scale <- seq(floor(min(mat, na.rm = TRUE)), 
                 round(max(mat, na.rm = TRUE), digits = 2) + 0.1, 
                 length.out = 6)
col_fun <- circlize::colorRamp2(mat_scale, 
                                pals::coolwarm(n = length(mat_scale)),
                                space = "LAB")

ht <- Heatmap(mat, col = col_fun, name = "JASMINE Score", 
              height = unit(11, "cm"), width = unit(5.5, "cm"), 
              column_names_gp = gpar(fontsize = 8), 
              row_names_gp = gpar(fontsize = 6),
              column_title =  "Lymphoid Cells scRNA-seq", 
              column_title_gp = gpar(fontsize = 12), 
              column_title_side = "top",
              row_dend_width = unit(5, "mm"), column_dend_height = unit(5, "mm"),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8), 
                                          labels_gp = gpar(fontsize = 6), 
                                          grid_height = unit(1, "cm"), 
                                          grid_width = unit(0.25, "cm")))

fid <- "Fig5F.lymphoid.JASMINE.heatmap.png"
png(fid, height = 13.5, width = 15, units = "cm", res = 200)
draw(ht)
invisible(dev.off())

# END