# reference: https://www.rhondabacher.com/scaffold-vignette.pdf

# pull data from La Manno et al (2016) -- https://doi.org/10.1016/j.cell.2016.09.027
fetch_lamanno_brain_data <- function() {
  sce <- scRNAseq::LaMannoBrainData(which = "human-embryo")
  sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
  sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with {scaffold}
  return(sce)
}

# pull data from Baron et al (2016) -- https://doi.org/10.1016/j.cels.2016.08.011
fetch_baron_pancreas_data <- function() {
  sce <- scRNAseq::BaronPancreasData(which = "human")
  sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
  sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with {scaffold}
  return(sce)
}

# pull data from Bastidas-Ponce et al (2019) -- https://doi.org/10.1242/dev.173849
fetch_bastidas_ponce_pancreas_data <- function(conda_env_name = "HPG_venv", conda_bin_path = "/apps/conda/22.11.1/condabin/conda") {
  # check inputs 
  if (!file.exists(conda_bin_path)) { stop("The provided conda binary doesn't exist.") }
  # fetch pancreatic endocrinogenesis data from `scVelo` Python package
  reticulate::use_condaenv(condaenv = conda_env_name,
                           conda = conda_bin_path, 
                           required = TRUE)
  scvelo <- reticulate::import("scvelo")
  adata <- scvelo$datasets$pancreas()
  sce <- zellkonverter::AnnData2SCE(adata = adata)
  sce@assays@data$X <- NULL
  sce@assays@data$counts <- sce@assays@data$spliced
  sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
  sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with {scaffold}
  return(sce)
}

# simulate data for one subject 
simulate_single_subject <- function(ref.dataset = NULL,
                                    perc.dyn.genes = NULL,
                                    n.cells = NULL, 
                                    n.knots = 2, 
                                    spline.degree = 2) {
  # check inputs
  if (is.null(ref.dataset) || is.null(perc.dyn.genes) || is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  
  # set up simulation parameters
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes <- ceiling((perc.dyn.genes / 0.8) * nrow(ref.dataset))  # make dynamic genes an 80% sample of the total pool of possible dynamic genes
  gene_means <- rowMeans(SingleCellExperiment::counts(ref.dataset))
  Q80 <- unname(stats::quantile(gene_means, 0.8))
  high_exp_genes <- rownames(ref.dataset)[gene_means > Q80]
  if (length(high_exp_genes) < n_possible_dyn_genes) {
    Q70 <- unname(stats::quantile(gene_means, 0.7))
    high_exp_genes <- rownames(ref.dataset)[gene_means > Q70]
    if (length(high_exp_genes) < n_possible_dyn_genes) {
      stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
    }
  }
  # make sure pool of possible dynamic genes has high expression
  possible_dyn_genes <- sample(high_exp_genes,
                               size = n_possible_dyn_genes, 
                               replace = FALSE)
  dyn_genes <- sample(possible_dyn_genes, 
                      n_dyn_genes, 
                      replace = FALSE)
  # set up B-spline knots 
  my_knots <- purrr::map(seq(n_dyn_genes), function(k) {
    knot_1 <- stats::runif(1, 0.1, 0.5)
    knot_2 <- stats::runif(1, 0.5, 0.9)
    while (abs(knot_1 - knot_2) < 0.1) {
      knot_2 <- stats::runif(1, 0.5, 0.9)
    }
    other_knots <- c()
    if (n.knots > 2) {
      other_knots <- stats::runif(n.knots - 2, 0.1, 0.9)
    }
    knot_mat <- matrix(c(knot_1, knot_2, other_knots), 
                       nrow = 1, 
                       ncol = n.knots)
    return(knot_mat)
  })
  my_knots <- purrr::reduce(my_knots, rbind)
  rownames(my_knots) <- dyn_genes
  
  # set up B-spline theta (between-knot directions)
  ncol_theta <- n.knots + spline.degree + 1
  theta_pop <- stats::rnorm(5 * n_dyn_genes * ncol_theta, 5, 5)
  my_theta <- matrix(sample(theta_pop, ncol_theta * n_dyn_genes), 
                     ncol = ncol_theta, 
                     nrow = n_dyn_genes, 
                     byrow = TRUE)
  rownames(my_theta) <- dyn_genes
  dynamic_params <- list(propGenes = perc.dyn.genes,
                         dynGenes = dyn_genes, 
                         degree = spline.degree,
                         knots = my_knots,
                         theta = my_theta)
  
  # simulate 10X dataset
  scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.dataset,
                                                          sceUMI = TRUE,
                                                          useUMI = TRUE,
                                                          protocol = "droplet",
                                                          numCells = n.cells,
                                                          popHet = c(1, 1),
                                                          useDynamic = dynamic_params)
  sim_data <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
  sim_data@assays@data$counts <- Matrix::Matrix(sim_data@assays@data$umi_counts, sparse = TRUE)
  sim_data@assays@data$umi_counts <- NULL
  sim_data <- sim_data[rownames(ref.dataset), ]  # original ordering
  
  # typical scran + scater pre-processing pipeline
  sim_data <- sim_data[rowSums(SingleCellExperiment::counts(sim_data) > 0) >= 3, colSums(SingleCellExperiment::counts(sim_data)) > 0] 
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data,
                             subset_row = top2k_hvgs, 
                             ncomponents = 50)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:30]
  sim_data <- scater::runUMAP(sim_data, 
                              dimred = "PCAsub", 
                              ncomponents = 2)
  g <- scran::buildSNNGraph(sim_data, 
                            use.dimred = "PCAsub", 
                            k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  SingleCellExperiment::colLabels(sim_data) <- factor(clusters)
  SummarizedExperiment::colData(sim_data) <- SummarizedExperiment::colData(sim_data) %>%
                                             as.data.frame() %>%
                                             dplyr::mutate(cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                                                           cell_time_normed = cell_time / max(cell_time)) %>%
                                             S4Vectors::DataFrame()
  return(sim_data)
}

# simulate data for multiple subjects
simulate_multi_subject <- function(ref.dataset = NULL,
                                   perc.dyn.genes = NULL,
                                   n.cells = NULL,
                                   perc.allocation = NULL,
                                   n.subjects = 6, 
                                   gene.dyn.threshold = 4, 
                                   n.knots = 2, 
                                   spline.degree = 2, 
                                   perc.overlap = 0.8) {
  # check inputs
  if (is.null(ref.dataset) || is.null(perc.dyn.genes) || is.null(n.cells)) { stop("You're missing vital parameters for simulate_scaffold().") }
  if (perc.dyn.genes <= 0) { stop("% dynamic genes need to be greater than zero.") }
  if (n.cells <= 0) { stop("Number of cells needs to be greater than zero.") }
  if (is.null(perc.allocation)) { stop("% allocation must be non-NULL.") }
  if (length(perc.allocation) != n.subjects) { stop("Each subject must have a % sample allocation value.") }
  
  # set up simulation parameters -- common across subjects
  obj_list <- vector("list", length = n.subjects)
  n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
  n_possible_dyn_genes <- ceiling((perc.dyn.genes / perc.overlap) * nrow(ref.dataset))  # make dynamic genes overlap by perc.overlap % between subject -- higher values -> more overlap
  gene_means <- rowMeans(SingleCellExperiment::counts(ref.dataset))
  Q80 <- unname(stats::quantile(gene_means, 0.8))
  high_exp_genes <- rownames(ref.dataset)[gene_means > Q80]
  if (length(high_exp_genes) < n_possible_dyn_genes) {
    Q70 <- unname(stats::quantile(gene_means, 0.7))
    high_exp_genes <- rownames(ref.dataset)[gene_means > Q70]
    if (length(high_exp_genes) < n_possible_dyn_genes) {
      stop("Your dataset has too few highly expressed genes to support the number of dynamic genes you want. Please reduce the % dynamic genes parameter.")
    }
  }
  # make sure pool of possible dynamic genes has high expression
  possible_dyn_genes <- sample(high_exp_genes,
                               size = n_possible_dyn_genes, 
                               replace = FALSE)
  # set up B-spline knots 
  my_knots <- purrr::map(seq(n_possible_dyn_genes), function(k) {
    knot_1 <- stats::runif(1, 0.1, 0.5)
    knot_2 <- stats::runif(1, 0.5, 0.9)
    while (abs(knot_1 - knot_2) < 0.1) {
      knot_2 <- stats::runif(1, 0.5, 0.9)
    }
    other_knots <- c()
    if (n.knots > 2) {
      other_knots <- stats::runif(n.knots - 2, 0.1, 0.9)
    }
    knot_mat <- matrix(c(knot_1, knot_2, other_knots), 
                       nrow = 1, 
                       ncol = n.knots)
    return(knot_mat)
  })
  my_knots <- purrr::reduce(my_knots, rbind)
  rownames(my_knots) <- possible_dyn_genes
  
  # set up B-spline theta (between-knot directions)
  ncol_theta <- n.knots + spline.degree + 1
  theta_pop <- stats::rnorm(5 * n_possible_dyn_genes * ncol_theta, 5, 5)
  my_theta <- matrix(sample(theta_pop, ncol_theta * n_possible_dyn_genes), 
                     ncol = ncol_theta, 
                     nrow = n_possible_dyn_genes, 
                     byrow = TRUE)
  rownames(my_theta) <- possible_dyn_genes
  
  # simulate 10X datasets
  for (s in seq(n.subjects)) {
    subject_n_cells <- ceiling(n.cells * perc.allocation[s])
    dyn_genes <- sample(possible_dyn_genes, 
                        n_dyn_genes, 
                        replace = FALSE)
    # fetch knots & theta values for each chosen dynamic gene from the overall pool of possible dynamic gene trends
    dyn_gene_knots <- my_knots[dyn_genes, ]
    dyn_gene_theta <- my_theta[dyn_genes, ]
    dynamic_params <- list(propGenes = perc.dyn.genes,
                           dynGenes = dyn_genes, 
                           degree = spline.degree,
                           knots = dyn_gene_knots,
                           theta = dyn_gene_theta)
    scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.dataset,
                                                            sceUMI = TRUE,
                                                            useUMI = TRUE,
                                                            protocol = "droplet",
                                                            numCells = subject_n_cells,
                                                            numGenes = nrow(ref.dataset), 
                                                            popHet = c(1, 1),
                                                            useDynamic = dynamic_params)
    sim_data <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
    sim_data@assays@data$counts <- Matrix::Matrix(sim_data@assays@data$umi_counts, sparse = TRUE)
    sim_data@assays@data$umi_counts <- NULL
    sim_data <- sim_data[rownames(ref.dataset), ]  # original ordering
    obj_list[[s]] <- sim_data
  }
  # combine & clean counts / metadata
  counts_mat <- purrr::map(obj_list, SingleCellExperiment::counts) %>%
                purrr::reduce(cbind)
  col_names <- purrr::map(seq(obj_list), \(i) paste0("P", i, "_", colnames(obj_list[[i]]))) %>% 
               purrr::reduce(c)
  row_names <- rownames(obj_list[[1]])  # common to all sims
  rownames(counts_mat) <- row_names
  colnames(counts_mat) <- col_names
  row_data <- purrr::map(obj_list, SummarizedExperiment::rowData) %>%
              purrr::reduce(cbind) %>%
              as.data.frame()
  colnames(row_data) <- paste0("geneStatus_P", 1:n.subjects)
  row_data <- dplyr::mutate(row_data, 
                            geneDynamic_n = rowSums(dplyr::across(dplyr::contains("geneStatus_"), \(x) x == "Dynamic")), 
                            geneStatus_overall = dplyr::if_else(geneDynamic_n >= gene.dyn.threshold, "Dynamic", "NotDynamic")) 
  row_data <- S4Vectors::DataFrame(row_data)
  col_data <- purrr::map(obj_list, SummarizedExperiment::colData) %>%
              purrr::reduce(rbind) %>%
              as.data.frame()
  subj_names <- purrr::map(seq(obj_list), \(i) rep(paste0("P", i), ncol(obj_list[[i]]))) %>% 
                purrr::reduce(c)
  cell_time_normed <- purrr::map(seq(obj_list), \(i) c(1:ncol(obj_list[[i]])) / ncol(obj_list[[i]])) %>% 
                      purrr::reduce(c)
  col_data <- dplyr::mutate(col_data,
                            subject = subj_names,
                            cell_time_normed = cell_time_normed)
  rownames(col_data) <- col_names
  col_data <- S4Vectors::DataFrame(col_data)
  sim_data <- SingleCellExperiment::SingleCellExperiment(list(counts = counts_mat))
  SummarizedExperiment::colData(sim_data) <- col_data
  SummarizedExperiment::rowData(sim_data) <- row_data
  sim_data <- sim_data[rownames(ref.dataset), ]  # original ordering -- just in case
  
  # process data w/ typical pipeline
  sim_data <- sim_data[rowSums(SingleCellExperiment::counts(sim_data) > 0) >= 3, colSums(SingleCellExperiment::counts(sim_data)) > 0] 
  sim_data <- scater::logNormCounts(sim_data)
  var_decomp <- scran::modelGeneVar(sim_data)
  top2k_hvgs <- scran::getTopHVGs(var_decomp, n = 2000)
  sim_data <- scater::runPCA(sim_data,
                             subset_row = top2k_hvgs, 
                             ncomponents = 50)
  SingleCellExperiment::reducedDim(sim_data, "PCAsub") <- SingleCellExperiment::reducedDim(sim_data, "PCA")[, 1:30]
  sim_data <- scater::runUMAP(sim_data, 
                              dimred = "PCAsub", 
                              ncomponents = 2)
  g <- scran::buildSNNGraph(sim_data, 
                            use.dimred = "PCAsub", 
                            k = 30)
  clusters <- igraph::cluster_louvain(graph = g)$membership
  SingleCellExperiment::colLabels(sim_data) <- factor(clusters)
  return(sim_data)
}
