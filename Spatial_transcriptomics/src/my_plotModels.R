theme_scLANE <- function(base.size = 12,
                         base.lwd = 0.75,
                         base.family = "sans") {
  ggplot2::theme_classic(base_size = base.size,
                         base_family = base.family,
                         base_line_size = base.lwd,
                         base_rect_size = base.lwd) +
    ggplot2::theme(strip.clip = "off",
                   strip.background = ggplot2::element_rect(linewidth = base.lwd))
}


my_plotModels <- function(test.dyn.res = NULL,
                       gene = NULL,
                       pt = NULL,
                       expr.mat = NULL,
                       size.factor.offset = NULL,
                       is.gee = FALSE,
                       is.glmm = FALSE,
                       id.vec = NULL,
                       cor.structure = "ar1",
                       ci.alpha = 0.05,
                       plot.null = TRUE,
                       plot.glm = TRUE,
                       plot.gam = TRUE,
                       plot.scLANE = TRUE,
                       filter.lineage = NULL,
                       gg.theme = theme_scLANE()) {
  # check inputs
  if (is.null(expr.mat) || is.null(pt) || is.null(gene) || is.null(test.dyn.res)) { stop("You forgot one or more of the arguments to plotModels().") }
  # get raw counts from SingleCellExperiment or Seurat object & transpose to cell x gene dense matrix
  if (inherits(expr.mat, "SingleCellExperiment")) {
    expr.mat <- as.matrix(BiocGenerics::counts(expr.mat))
  } else if (inherits(expr.mat, "Seurat")) {
    expr.mat <- as.matrix(Seurat::GetAssayData(expr.mat,
                                               slot = "counts",
                                               assay = Seurat::DefaultAssay(expr.mat)))
  } else if (inherits(expr.mat, "dgCMatrix")) {
    expr.mat <- as.matrix(expr.mat)
  }
  if (!(inherits(expr.mat, "matrix") || inherits(expr.mat, "array"))) { stop("Input expr.mat must be coerceable to a matrix of integer counts.") }
  expr.mat <- t(expr.mat)
  # generate parameters for CIs
  Z <- stats::qnorm(ci.alpha / 2, lower.tail = FALSE)
  # select sublist for gene of interest
  td_res <- test.dyn.res[[gene]]
  # make sure lineages are named nicely
  colnames(pt) <- paste0("Lineage_", LETTERS[1:ncol(pt)])
  # create base list w/ elements being lineage-specific dataframes
  counts_df_list <- purrr::map2(pt,
                                LETTERS[1:ncol(pt)],
                                \(x, y) {
                                  mod_df <- data.frame(CELL = rownames(pt)[!is.na(x)],
                                                       LINEAGE = y,
                                                       PT = x[!is.na(x)],
                                                       COUNT = expr.mat[!is.na(x), gene])
                                  if (is.gee || is.glmm) {
                                    mod_df <- dplyr::mutate(mod_df, ID = id.vec[!is.na(x)])
                                  } else {
                                    mod_df <- dplyr::mutate(mod_df, ID = NA_character_)
                                  }
                                  if (!is.null(size.factor.offset)) {
                                    mod_df <- dplyr::mutate(mod_df, CELL_OFFSET = unname(size.factor.offset[!is.na(x)]))
                                  }
                                  return(mod_df)
                                })
  # create list of dataframes w/ predicted values, standard errors, CIs for null, GLM, GAM, & MARGE models
  counts_df_list <- purrr::map(counts_df_list, \(x) dplyr::relocate(x, ID, CELL, LINEAGE, COUNT, PT)) %>%
                    purrr::map2(.y = td_res, \(.x, .y) dplyr::mutate(.x,
                                                                     RESP_MARGE = .y$MARGE_Preds$marge_link_fit,
                                                                     SE_MARGE = .y$MARGE_Preds$marge_link_se,
                                                                     RESP_NULL = .y$Null_Preds$null_link_fit,
                                                                     SE_NULL = .y$Null_Preds$null_link_se)) %>%
                    purrr::map(function(x) {
                      if (is.gee) {
                        theta_hat <- MASS::theta.mm(y = x$COUNT,
                                                    mu = mean(x$COUNT),
                                                    dfr = nrow(x) - 1)
                        gee_formula <- "COUNT ~ PT"
                        if (!is.null(size.factor.offset)) {
                          gee_formula <- paste0(gee_formula, " + offset(log(1 / CELL_OFFSET))")
                        }
                        gee_formula <- stats::as.formula(gee_formula)
                        glm_mod <- geeM::geem(gee_formula,
                                              id = ID,
                                              data = x,
                                              family = MASS::negative.binomial(theta_hat, link = log),
                                              corstr = cor.structure,
                                              scale.fix = FALSE,
                                              sandwich = TRUE)
                        robust_vcov_mat <- as.matrix(glm_mod$var)
                        glm_preds <- data.frame(fit = predict(glm_mod),
                                                se.fit = sqrt(apply((tcrossprod(glm_mod$X, robust_vcov_mat)) * glm_mod$X, 1, sum)))
                      } else if (is.glmm) {
                        if (is.null(size.factor.offset)) {
                          glm_mod <- glmmTMB::glmmTMB(COUNT ~ PT + (1 | ID) + (1 + PT | ID),
                                                      data = x,
                                                      family = glmmTMB::nbinom2(link = "log"),
                                                      se = TRUE,
                                                      REML = FALSE)
                        } else {
                          glm_mod <- glmmTMB::glmmTMB(COUNT ~ PT + (1 | ID) + (1 + PT | ID),
                                                      data = x,
                                                      family = glmmTMB::nbinom2(link = "log"),
                                                      offset = log(1 / x$CELL_OFFSET),
                                                      se = TRUE,
                                                      REML = FALSE)
                        }
                        glm_preds <- data.frame(predict(glm_mod, type = "link", se.fit = TRUE)[1:2])
                      } else {
                        glm_formula <- "COUNT ~ PT"
                        if (!is.null(size.factor.offset)) {
                          glm_formula <- paste0(glm_formula, " + offset(log(1 / CELL_OFFSET))")
                        }
                        glm_formula <- stats::as.formula(glm_formula)
                        glm_mod <- MASS::glm.nb(glm_formula,
                                                data = x,
                                                x = FALSE,
                                                y = FALSE,
                                                method = "glm.fit2",
                                                link = log,
                                                init.theta = 1)
                        glm_preds <- data.frame(stats::predict(glm_mod, type = "link", se.fit = TRUE)[1:2])
                      }
                      x %<>% dplyr::mutate(RESP_GLM = glm_preds$fit,
                                           SE_GLM = glm_preds$se.fit)
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      if (is.null(size.factor.offset)) {
                        if (is.glmm) {
                          gam_mod <- nbGAM(expr = x$COUNT,
                                           pt = x$PT,
                                           id.vec = x$ID)
                        } else {
                          gam_mod <- nbGAM(expr = x$COUNT, pt = x$PT)
                        }
                      } else {
                        if (is.glmm) {
                          gam_mod <- nbGAM(expr = x$COUNT,
                                           Y.offset = x$CELL_OFFSET,
                                           pt = x$PT,
                                           id.vec = x$ID)
                        } else {
                          gam_mod <- nbGAM(expr = x$COUNT,
                                           Y.offset = x$CELL_OFFSET,
                                           pt = x$PT)
                        }
                      }
                      gam_preds <- data.frame(predict(gam_mod, type = "link", se.fit = TRUE)[1:2])
                      x <- dplyr::mutate(x,
                                         RESP_GAM = gam_preds$fit,
                                         SE_GAM = gam_preds$se.fit)
                      return(x)
                    }) %>%
                    purrr::map(function(x) {
                      dplyr::select(x, CELL, CELL_OFFSET, ID, LINEAGE, COUNT, PT) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("RESP")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("RESP"),
                                                            values_to = "RESPONSE",
                                                            names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("RESP_", "", MODEL))),
                                       by = c("CELL" = "CELL",
                                              "ID" = "ID",
                                              "LINEAGE" = "LINEAGE")) %>%
                      dplyr::left_join((dplyr::select(x, CELL, ID, LINEAGE, dplyr::contains("SE")) %>%
                                        tidyr::pivot_longer(cols = dplyr::contains("SE"),
                                                            values_to = "SE",
                                                            names_to = "MODEL") %>%
                                        dplyr::mutate(MODEL = gsub("SE_", "", MODEL))),
                                       by = c("CELL" = "CELL",
                                              "ID" = "ID",
                                              "LINEAGE" = "LINEAGE",
                                              "MODEL" = "MODEL"))
    })
  counts_df <- purrr::reduce(counts_df_list, rbind) %>%
               dplyr::mutate(MODEL = dplyr::case_when(MODEL == "NULL" ~ "Intercept-only",
                                                      MODEL == "MARGE" ~ "scLANE",
                                                      TRUE ~ MODEL),
                             MODEL = factor(MODEL, levels = c("Intercept-only", "GLM", "GAM", "scLANE")),
                             PRED = exp(RESPONSE),
                             CI_LL = exp(RESPONSE - Z * SE),
                             CI_UL = exp(RESPONSE + Z * SE))
  if (!is.null(size.factor.offset)) {
    if (is.gee) {
      counts_df <- dplyr::rowwise(counts_df) %>%
                   dplyr::mutate(COUNT = COUNT * CELL_OFFSET,
                                 dplyr::across(c(PRED, CI_LL, CI_UL), \(x) {
                                   if (as.character(MODEL) == "GAM") {
                                     x * CELL_OFFSET
                                   } else {
                                     x
                                   }
                                 }))
    } else {
      counts_df <- dplyr::mutate(counts_df, dplyr::across(c(COUNT, PRED, CI_LL, CI_UL), \(x) x * CELL_OFFSET))
    }
  }
  # add conditional filters here
  if (!plot.null) {
    counts_df %<>% dplyr::filter(MODEL != "Intercept-only")
  }
  if (!plot.glm) {
    counts_df %<>% dplyr::filter(MODEL != "GLM")
  }
  if (!plot.gam) {
    counts_df %<>% dplyr::filter(MODEL != "GAM")
  }
  if (!plot.scLANE) {
    counts_df %<>% dplyr::filter(MODEL != "scLANE")
  }
  if (!is.null(filter.lineage)) {
    counts_df %<>% dplyr::filter(!LINEAGE %in% filter.lineage)
  }
  # change model labels as necessary
  if (is.gee) {
    counts_df %<>% dplyr::mutate(MODEL = dplyr::if_else(as.character(MODEL) == "GLM", "GEE", as.character(MODEL)),
                                 MODEL = factor(MODEL, levels = c("Intercept-only", "GEE", "GAM", "scLANE")))
  }
  if (is.glmm) {
    counts_df %<>% dplyr::mutate(MODEL = dplyr::if_else(as.character(MODEL) == "GLM", "GLMM", as.character(MODEL)),
                                 MODEL = factor(MODEL, levels = c("Intercept-only", "GLMM", "GAM", "scLANE")))
  }
  # generate plot
  if (is.glmm) {
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, group = ID)) +
         ggplot2::geom_point(size = 0.85, color = alpha("black", 0.5))
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_nested_wrap(~paste0("Lineage ", LINEAGE) + MODEL + ID,
                                        nrow = length(levels(counts_df$MODEL)),
                                        strip = ggh4x::strip_nested(clip = "off",
                                                                    background_x = list(ggplot2::element_rect(linewidth = gg.theme$line$linewidth))))
    } else {
      p <- p + ggplot2::facet_wrap(~paste0("Lineage ", LINEAGE) + MODEL + ID)
    }
    p <- p +
         ggplot2::geom_line(mapping = ggplot2::aes(y = PRED),
                            linewidth = 0.75,
                            color = "red2",
                            show.legend = FALSE) +
         ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = CI_LL, ymax = CI_UL),
                              fill = "grey50",
                              alpha = 0.4,
                              linewidth = 0,
                              show.legend = FALSE) +
         ggplot2::scale_y_continuous(labels = scales::label_comma()) +
         ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
         ggplot2::labs(x = "Pseudotime",
                       y = "Expression",
                       color = "Subject",
                       fill = "Subject",
                       title = gene) +
         gg.theme +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "italic"),
                        strip.clip = "off",
                        strip.background = ggplot2::element_rect(linewidth = gg.theme$line$linewidth)) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)))
  } else {
    p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = PT, y = COUNT, color = LINEAGE)) +
      ggplot2::geom_point(size = 0.85, color = alpha("black", 0.5),
                             show.legend = ifelse(ncol(pt) > 1, TRUE, FALSE))
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::facet_nested_wrap(~paste0("Lineage ", LINEAGE) + MODEL,
                                        strip = ggh4x::strip_nested(clip = "off",
                                                                    background_x = list(ggplot2::element_rect(linewidth = gg.theme$line$linewidth))))
    } else {
      p <- p + ggplot2::facet_wrap(~paste0("Lineage ", LINEAGE) + MODEL)
    }
    p <- p +
         ggplot2::geom_line(mapping = ggplot2::aes(y = PRED),
                            linewidth = 0.75,
                            color = "red2") +
         ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = CI_LL, ymax = CI_UL),
                              alpha = 0.4,
                              linewidth = 0,
                              color = "grey") +
         ggplot2::scale_y_continuous(labels = scales::label_comma()) +
         ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
         ggplot2::labs(x = "Pseudotime",
                       y = "Expression",
                       color = "Lineage",
                       fill = "Lineage",
                       title = gene) +
         gg.theme +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "italic"),
                        strip.clip = "off",
                        strip.background = ggplot2::element_rect(linewidth = gg.theme$line$linewidth)) +
         ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1)))
  }
  return(p)
}
