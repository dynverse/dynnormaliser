# lsei import is for scran

#' State-of-the-art preprocessing and normalisation
#'
#' Based on \url{https://f1000research.com/articles/5-2122/v2} and \url{https://www.bioconductor.org/help/workflows/simpleSingleCell/}.
#'
#' @param counts The counts matrix, with features in columns
#' @param filter_cells Whether the cells have to be filtered
#' @param filter_features Whether the features have to be filtered
#' @param filter_hvg Whether to filter on highly variable features
#' @param normalisation How to normalise
#' @param has_spike Does this contain spike-ins, for which the feature names are preseded by ERCC
#' @param verbose Whether to add plots
#' @param nmads Number of median deviations for filtering outlier cells
#' @param min_ave_expression Minimal average expression of a feature
#' @param hvg_fdr FDR feature filtering cutoff
#' @param hvg_bio Biological feature filtering cutoff
#' @param min_variable_fraction Minimal number of variable features to retain
#'
#' @importFrom Matrix t rowMeans
#' @importFrom scater calculateQCMetrics isOutlier normalise
#' @importFrom SingleCellExperiment isSpike SingleCellExperiment counts sizeFactors logcounts
#' @importFrom scran computeSumFactors computeSpikeFactors trendVar decomposeVar
#' @importFrom limSolve lsei
#' @importFrom stats sd
#' @importFrom aroma.light aroma.light
#' @export
normalise_filter_counts <- function(
  counts,
  filter_cells = TRUE,
  filter_features = TRUE,
  filter_hvg = TRUE,
  normalisation = "scran_size_factors",
  has_spike = any(grepl("^ERCC", colnames(counts))),
  verbose = FALSE,
  nmads = 3,
  min_ave_expression = 0.02,
  hvg_fdr = 0.05,
  hvg_bio = 0.5,
  min_variable_fraction = 0.1
) {
  if (verbose) {
    requireNamespace("grDevices")
    requireNamespace("ggplot2")
    requireNamespace("graphics")
    requireNamespace("KernSmooth")
  }

  normalisation_plots <- list()
  normalisation_steps <- tibble()

  ########################################
  # Create data object
  ########################################

  # convert to integer
  counts <- floor(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Matrix::t(counts)))

  # mitochondrial
  mitochondrial <- grepl("^(mt|MT|Mt)-", rownames(sce))
  has_mito <- any(mitochondrial)

  # spike ins
  spike <- grepl("^ERCC", rownames(sce))
  has_spike <- any(spike)

  # calculate qc metrics
  feature_controls <- list()

  if (has_mito) feature_controls$Mt <- mitochondrial
  if (has_spike) feature_controls$ERCC <- grepl("^ERCC", rownames(sce))

  sce <- scater::calculateQCMetrics(
    sce,
    feature_controls = feature_controls,
    compact = TRUE
  )

  if (has_spike) {
    SingleCellExperiment::isSpike(sce, "ERCC") <- spike
  }

  # plots
  if (verbose) {
    normalisation_steps <- tribble(
      ~type, ~nfeatures, ~ncells,
      "original", nrow(sce), ncol(sce)
    )
    print(glue::glue("Original: features - {nrow(sce)} Cells - {ncol(sce)}"))

    normalisation_plots$library <-
      ggplot(as.data.frame(sce$scater_qc$all)) +
      geom_histogram(aes(total_counts)) +
      scale_x_continuous(limits = c(0, NA))

    if (has_spike) {
      normalisation_plots$spike <-
        ggplot(as.data.frame(sce$scater_qc$feature_control_ERCC)) +
        geom_histogram(aes(pct_counts)) +
        scale_x_continuous(limits = c(0, 100))
    }

    if (has_mito) {
      normalisation_plots$mito <-
        ggplot(as.data.frame(sce$scater_qc$feature_control_Mt)) +
        geom_histogram(aes(pct_counts)) +
        scale_x_continuous(limits = c(0, 100))
    }

  }

  ########################################
  # Filter cells
  ########################################

  if (filter_cells) {
    total_counts <- sce$scater_qc$all$log10_total_counts
    total_features <- sce$scater_qc$all$log10_total_features_by_counts
    pct_counts_Mt <- sce$scater_qc$feature_control_Mt$pct_counts
    pct_counts_ERCC <- sce$scater_qc$feature_control_ERCC$pct_counts

    mito_drop <- rep(FALSE, length(total_counts))
    spike_drop <- rep(FALSE, length(total_counts))

    libsize_drop <- scater::isOutlier(total_counts, nmads = nmads, type = "lower", log = TRUE)
    feature_drop <- scater::isOutlier(total_features, nmads = nmads, type = "lower", log = TRUE)
    if (has_mito) mito_drop <- scater::isOutlier(pct_counts_Mt, nmads = nmads, type = "higher")
    if (has_spike) spike_drop <- scater::isOutlier(pct_counts_ERCC, nmads = nmads, type = "higher")

    if (verbose) {
      tibble(sum(mito_drop), sum(spike_drop), sum(libsize_drop), sum(feature_drop)) %>% print()
    }

    sce <- sce[,!(libsize_drop | feature_drop | mito_drop | spike_drop)]

    if (verbose) {
      normalisation_steps <- normalisation_steps %>%
        add_row(type = "cell_quality_filtering", nfeatures = nrow(sce), ncells = ncol(sce))
      print(glue::glue("Cell filter: features - {nrow(sce)} Cells - {ncol(sce)}"))
    }
  }

  ########################################
  # Filter features
  ########################################

  if (filter_features) {
    ave_counts <- Matrix::rowMeans(BiocGenerics::counts(sce))
    keep <- ave_counts >= min_ave_expression

    sce <- sce[keep,]

    if (verbose) {
      normalisation_plots$ave_counts <- tibble(ave_counts = ave_counts) %>%
        ggplot() +
        geom_histogram(aes(ave_counts)) +
        scale_x_log10() +
        geom_vline(xintercept = min_ave_expression)

      top_features <- ave_counts %>% sort() %>% tail(20) %>% names()
      counts_top_features <- SingleCellExperiment::counts(sce[top_features]) %>%
        reshape2::melt(varnames = c("feature", "cell"), value.name = "count") %>%
        dplyr::mutate(feature = factor(feature, levels = top_features))

      avecounts_top_features <- counts_top_features %>%
        group_by(feature) %>%
        summarise(count = mean(count))

      normalisation_plots$top_counts <- counts_top_features %>%
        ggplot(aes(feature, count + 1)) +
        geom_point(shape = "|") +
        geom_point(data = avecounts_top_features, color = "blue", size = 4) +
        scale_y_log10() +
        coord_flip()

      normalisation_steps <- normalisation_steps %>%
        add_row(type = "feature_expression_filtering", nfeatures = nrow(sce), ncells = ncol(sce))
      print(glue::glue("feature filter: features - {nrow(sce)} Cells - {ncol(sce)}"))
    }
  }

  ########################################
  # Normalise
  ########################################

  if (normalisation == "scran_size_factors") {
    # fix for small datasets with very low number of cells
    if (ncol(sce) >= 100) {
      sizes <- c(20, 40, 60, 80)
    } else {
      sizes <- round(ncol(sce) / 5)
    }

    # In some cases, the computation of these sum factors can lead to very strange size factors, which are not at all correlated with the total counts in each cell
    # If that is the case, we for now fall back to scater::librarySizeFactors
    sce <- scran::computeSumFactors(sce, sizes = sizes, min.mean = 0.1, positive = TRUE)
    sce <- sce[, SingleCellExperiment::sizeFactors(sce) > 0] # as mentioned in the scran documentation, ensure that size factors are higher than 0

    if(has_spike) {
      sce <- scran::computeSpikeFactors(sce, type = "ERCC", general.use = FALSE)
      if(any(is.na(sizeFactors(sce, type = "ERCC")))) {
        warning("Some cells do not have any spike-ins, this will cause an error further away. Remove spike-ins.")
      }
    }

    if (cor(sce$scater_qc$all$total_counts, sizeFactors(sce)) < 0.5) {
      warning("Low correlation between sumFactors and total_counts, falling back to scater::librarySizeFactors")
      SingleCellExperiment::sizeFactors(sce) <- scater::librarySizeFactors(sce)
    }

    sce <- scater::normalise(sce)

    if (verbose) {
      normalisation_plots$size_factors <- tibble(
        size_factor = sizeFactors(sce),
        library_size = sce$scater_qc$all$total_counts
      ) %>%
        ggplot(aes(size_factor, library_size)) +
        geom_point() +
        geom_smooth(method = "lm")
    }
  } else {
    stop("Normalisation not supported")
  }

  if (verbose) {
    normalisation_steps <- normalisation_steps %>%
      add_row(type = "normalisation", nfeatures = nrow(sce), ncells = ncol(sce))
    print(glue::glue("Normalised: features - {nrow(sce)} Cells - {ncol(sce)}"))
  }

  ########################################
  # Select highly variable features
  ########################################

  if (filter_hvg) {
    var_fit <- scran::trendVar(sce, method = "spline", use.spikes = has_spike) # requires aroma.light
    var_out <- scran::decomposeVar(sce, var_fit) %>% as.data.frame()

    if (verbose) {
      if (has_spike) {
        var_out$spike <- SingleCellExperiment::isSpike(sce)
      } else {
        var_out$spike <- FALSE
      }

      normalisation_plots$variance_mean <- var_out %>%
        ggplot(aes(mean, total)) +
          geom_point(aes(color = spike)) +
          geom_smooth()

      normalisation_plots$fdr_bio <- var_out %>%
        ggplot(aes(FDR, bio)) +
          geom_point() +
          geom_hline(yintercept = hvg_bio, color = "red") +
          geom_vline(xintercept = hvg_fdr, color = "red")
    }

    var_out <- var_out[order(var_out$bio, decreasing = TRUE),]
    hvg_out <- var_out[which(var_out$FDR <= hvg_fdr & var_out$bio >= hvg_bio),]
    if(nrow(hvg_out) < min_variable_fraction * nrow(var_out)) {
      hvg_out <- var_out[seq(1, ceiling(min_variable_fraction * nrow(var_out))), ]
    }

    sce <- sce[rownames(hvg_out),]

    if (verbose) {
      normalisation_steps <- normalisation_steps %>%
        add_row(type = "feature_variability_filtering", nfeatures = nrow(sce), ncells = ncol(sce))
      print(glue::glue("Variable features filtered: features - {nrow(sce)} Cells - {ncol(sce)}"))
    }
  }

  expr_norm_filt <- SingleCellExperiment::logcounts(sce) %>% Matrix::t()
  counts_filt <- counts[rownames(expr_norm_filt),colnames(expr_norm_filt)]

  ########################################
  # Iterative filtering on variability
  ########################################
  repeat {
    feature_sds <- counts_filt %>% apply(2, stats::sd)
    cell_sds <- counts_filt %>% apply(1, stats::sd)

    features_filtered <- which(feature_sds > 0, useNames = TRUE)
    cells_filtered <- which(cell_sds > 0, useNames = TRUE)
    expr_norm_filt <- expr_norm_filt[cells_filtered, features_filtered]
    counts_filt <- counts_filt[cells_filtered, features_filtered]

    if((min(c(feature_sds, 1), na.rm = TRUE) > 0) && (min(c(cell_sds, 1), na.rm = TRUE) > 0)){
      break
    }
  }

  if (verbose) {
    normalisation_steps <- normalisation_steps %>%
      add_row(type = "final_filtering", nfeatures = ncol(expr_norm_filt), ncells = nrow(expr_norm_filt))
    print(glue::glue("Final filtering: features - {ncol(expr_norm_filt)} Cells - {nrow(expr_norm_filt)}"))
  }

  ########################################
  # Output
  ########################################

  if(verbose) {
    type <- NULL # satisfy r cmd check

    normalisation_plots$n_retained <-
      normalisation_steps %>%
      mutate(type = factor(type, levels = rev(type))) %>%
      gather("dimension", "n", -type) %>%
      ggplot2::ggplot() +
      ggplot2::geom_bar(ggplot2::aes_string("type", "n", fill = "dimension"), position = "dodge", stat = "identity") +
      ggplot2::facet_wrap(~dimension, scales = "free_x") +
      ggplot2::coord_flip()
  } else {
    normalisation_steps <- NULL
  }

  lst(
    expression = expr_norm_filt,
    counts = counts_filt,
    normalisation_plots,
    info = lst(
      has_spike,
      has_mito,
      normalisation_steps
    )
  )
}
