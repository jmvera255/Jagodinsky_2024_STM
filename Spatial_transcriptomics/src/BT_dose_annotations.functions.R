
# BT_dose_annotations.functions.R

circle_test <- function(x, y, pt_source, r) {
  h <- pt_source$pxl_col_in_fullres
  k <- pt_source$pxl_row_in_fullres
  
  r^2 >= (x - h)^2 + (y-k)^2
}

calc_distance <- function(x, y, pt_source) {
  h <- pt_source$pxl_col_in_fullres
  k <- pt_source$pxl_row_in_fullres
  
  r <- sqrt((x - h)^2 + (y-k)^2)
  
  return(r)
}

get_pt_source <- function(spatial_df, spe_coldata, pt_source) {
  spe_coldata <- left_join(spe_coldata, spatial_df, by = "barcode")
  row_barcode <- filter(pt_source, pt_source == "row") %>% pull(Barcode)
  row <- filter(spe_coldata, barcode %in% row_barcode) %>% pull(array_row)
  
  col_barcode <- filter(pt_source, pt_source == "col") %>% pull(Barcode)
  col <- filter(spe_coldata, barcode == col_barcode[1]) %>% pull(array_col)
  
  pt_source <- rbind(filter(spe_coldata, array_col == col & array_row == row[1]), 
                     filter(spe_coldata, array_col == col & array_row == row[2])) %>%
    select(pxl_col_in_fullres, pxl_row_in_fullres, barcode)
  
  return(pt_source)
}

define_regions <- function(dist, px_per_um, spatial_df, pt_source){
  pxl_dist <- dist*px_per_um
  min_row <- pt_source$pxl_row_in_fullres-pxl_dist
  max_row <- pt_source$pxl_row_in_fullres+pxl_dist
  min_col <- pt_source$pxl_col_in_fullres-pxl_dist
  max_col <- pt_source$pxl_col_in_fullres+pxl_dist
  test_spots <- filter(spatial_df, pxl_col_in_fullres < max_col & pxl_col_in_fullres > min_col)
  test_spots <- filter(test_spots, pxl_row_in_fullres < max_row & pxl_row_in_fullres > min_row)
  test_spots <- mutate(test_spots, 
                       test = unlist(purrr::map2(pxl_col_in_fullres, pxl_row_in_fullres, 
                                                 function(x,y) circle_test(x, y, 
                                                                           pt_source,
                                                                           pxl_dist))))
  test_spots <- left_join(spatial_df, select(test_spots, barcode, test), by = "barcode")
  test_spots$test[is.na(test_spots$test)] <- FALSE
  return(test_spots)
}