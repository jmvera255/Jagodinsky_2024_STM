
# BayesSpace.functions.R

prep_spe <- function(spe) {
  # replace raw counts with SpotClean decontaminated counts
  assays(spe)$counts <- assays(spe)$decont
  keep <- which(colData(spe)$spot_QC == "keep")
  # remove poor QC spots
  spe <- spe[,keep]
  
  # add appropriately labeled colData so I don't have to revise
  # the BayesSpace code
  i <- grep("array_row", names(colData(spe)))
  names(colData(spe))[i] <- "row"
  i <- grep("array_col", names(colData(spe)))
  names(colData(spe))[i] <- "col"
  #colData(spe)$row <- colData(spe)$array_row
  #colData(spe)$col <- colData(spe)$array_col
  
  set.seed(1234)
  spe <- suppressMessages(BayesSpace::spatialPreprocess(spe, platform = "Visium", n.PCs = 15, log.normalize = TRUE))
  set.seed(1234)
  spe <- suppressMessages(BayesSpace::qTune(spe, qs=2:18, platform = "Visium", d = 15))
  
  return(spe)
}

estimate_q <- function(spe) {
  sample_id <- unique(spe$sample_id)
  
  p <- qPlot(spe) + 
    scale_x_continuous(name = "# of clusters", 
                       breaks = 2:18, 
                       labels = as.character(2:18)) +
    ggtitle(sample_id)
  return(p)
}

get_bayesspace_clusters <- function(spe, q) {
  spe <- suppressMessages(spatialCluster(spe, q = q, platform = "Visium", d = 15, 
                                                    init.method = "mclust", model = "t", 
                                                    gamma = 2, nrep = 50000, burn.in = 1000, 
                                                    #gamma = 2, nrep = 2000, burn.in = 100, 
                                                    save.chain = FALSE))
  clusters <- data.frame("BayesSpace_clusters" = spe$spatial.cluster)
  rownames(clusters) <- colnames(spe)
  return(clusters)
}





# END