
# Seurat_transferData.functions.R
convert_to_seurat_spatial <- function(spe) {
  sample_id <- unique(spe$sample_id)
  image.dir <- paste0(sample_id[1], "/outs/spatial")
  
  idx <- which(spe$spot_QC == "keep")
  spe <- spe[,idx]
  spot_metadata <- colData(spe)
  spot_metadata$barcode <- colnames(spe)
  rownames(spot_metadata) <- colnames(spe)
  
  decont <- assays(spe)$decont

  rownames(decont) <- spe@rowRanges@elementMetadata@listData$symbol
  
  image <- png::readPNG(source = file.path(image.dir,
                                           image.name = "tissue_lowres_image.png"))
  
  scale.factors <- jsonlite::fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  
  tissue.positions <- cbind.data.frame(colData(spe), spatialCoords(spe))
  tissue.positions <- dplyr::select(tissue.positions, in_tissue, 
                                    array_row, array_col, 
                                    pxl_row_in_fullres, pxl_col_in_fullres)
  colnames(tissue.positions) <- c('tissue', 'row', 'col', 
                                  'imagerow', 'imagecol')
  rownames(tissue.positions) <- rownames(colData(spe))
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  
  visium <- new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef),
    coordinates = tissue.positions,
    spot.radius = spot.radius)
  
  #meta <- cbind.data.frame(colData(spe)) %>% select(matches("spatial"))
  
  object <- CreateSeuratObject(counts = decont, assay = "Spatial")
  visium <- visium[Cells(x = object)]
  DefaultAssay(object = visium) <- "Spatial"
  object[["slice1"]] <- visium
  object$orig.ident <- unique(spe$sample_id)
  object[["Spatial"]][["ensembl_gene_id"]] <- rownames(spe)
  object@meta.data <- spot_metadata
  
  
  
  return(object)
}
