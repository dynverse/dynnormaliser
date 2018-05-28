# lsei import is for scran

#' State-of-the-art preprocessing and normalisation
#'
#' Based on \url{https://f1000research.com/articles/5-2122/v2} and \url{https://www.bioconductor.org/help/workflows/simpleSingleCell/}.
#'
#' @param counts The counts matrix, with genes in columns
#' @param has_spike Does this contain spike-ins, for which the gene names are preseded by ERCC
#' @param verbose Whether to add plots
#' @param nmads Number of median deviations for filtering outlier cells
#' @param expressed_in_n_cells Percentage of minimal number of cells a gene has to be expressed
#' @param min_ave_expression Minimal average expression of a gene
#' @param filter_hvg Whether to filter out highly variable genes
#' @param hvg_fdr FDR gene filtering cutoff
#' @param hvg_bio Biological gene filtering cutoff
#' @param min_variable_fraction Minimal number of variable genes to retain
#'
#' @importFrom Matrix t rowMeans
#' @importFrom scater calculateQCMetrics isOutlier normalise plotExplanatoryVariables plotExpression plotQC nexprs
#' @importFrom SingleCellExperiment isSpike SingleCellExperiment
#' @importFrom BiocGenerics counts sizeFactors
#' @importFrom Biobase pData
#' @importFrom scran computeSumFactors computeSpikeFactors trendVar decomposeVar
#' @importFrom limSolve lsei
#' @importFrom stats sd
#' @importFrom aroma.light aroma.light
#' @export
normalise_filter_counts <- function(
  counts,
  has_spike = any(grepl("^ERCC", colnames(counts))),
  verbose = FALSE,
  nmads = 3,
  expressed_in_n_cells = 0.05,
  min_ave_expression = 0.05,
  filter_hvg = TRUE,
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

  ########################################
  # Create data object
  ########################################

  counts <- round(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=Matrix::t(counts)))

  mitochondrial <- grepl("^(mt|MT|Mt)-", rownames(sce))
  has_mito <- any(mitochondrial)
  feature_controls <- list()
  if (has_mito) feature_controls$Mt <- mitochondrial
  if (has_spike) feature_controls$ERCC <- grepl("^ERCC", rownames(sce))

  sce <- scater::calculateQCMetrics(sce, feature_controls = feature_controls)

  if (has_spike) {
    is_spike <- grepl("^ERCC", rownames(sce))
    SingleCellExperiment::isSpike(sce, "ERCC") <- is_spike
    summary(is_spike)
  }

  if (verbose) {
    normalisation_steps <- tribble(
      ~type, ~ngenes, ~ncells,
      "original", nrow(sce), ncol(sce)
    )
    print(glue::glue("Original: Genes - {nrow(sce)} Cells - {ncol(sce)}"))
  }

  if (verbose) {
    graphics::par(mfrow=c(1,2))
    graphics::hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",
                   breaks=20, col="grey80", ylab="Number of cells")
    graphics::hist(sce$total_features, xlab="Number of expressed genes", main="",
                   breaks=20, col="grey80", ylab="Number of cells")
    graphics::par(mfrow=c(1, 1))
    normalisation_plots$library <- grDevices::recordPlot()
  }

  if (verbose) {
    graphics::par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
    graphics::hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",
                   breaks=20, col="grey80", ylab="Number of cells")
    graphics::hist(sce$total_features, xlab="Number of expressed genes", main="",
                   breaks=20, col="grey80", ylab="Number of cells")
    if (has_spike) graphics::hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)",
                                  ylab="Number of cells", breaks=20, main="", col="grey80")
    if (has_mito) graphics::hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)",
                                 ylab="Number of cells", breaks=20, main="", col="grey80")
    graphics::par(mfrow=c(1, 1))
    normalisation_plots$cell_quality <- grDevices::recordPlot()
  }

  ########################################
  # Filter cells
  ########################################
  mito_drop <- rep(FALSE, length(sce$total_counts))
  spike_drop <- rep(FALSE, length(sce$total_counts))

  libsize_drop <- scater::isOutlier(sce$total_counts, nmads=nmads, type="lower", log=TRUE)
  feature_drop <- scater::isOutlier(sce$total_features, nmads=nmads, type="lower", log=TRUE)
  if (has_mito) mito_drop <- scater::isOutlier(sce$pct_counts_Mt, nmads=nmads, type="higher")
  if (has_spike) spike_drop <- scater::isOutlier(sce$pct_counts_ERCC, nmads=nmads, type="higher")

  if (verbose) {
    tibble(sum(mito_drop), sum(spike_drop), sum(libsize_drop), sum(feature_drop)) %>% print()
  }

  sce <- sce[,!(libsize_drop | feature_drop | mito_drop | spike_drop)]

  if (verbose) {
    normalisation_steps <- normalisation_steps %>%
      add_row(type = "cell_quality_filtering", ngenes = nrow(sce), ncells = ncol(sce))
    print(glue::glue("Cell filter: Genes - {nrow(sce)} Cells - {ncol(sce)}"))
  }

  ########################################
  # Filter genes
  ########################################

  ave_counts <- Matrix::rowMeans(BiocGenerics::counts(sce))
  keep <- ave_counts >= min_ave_expression

  if (verbose) {
    fontsize <- ggplot2::theme(
      axis.text = ggplot2::element_text(size=12),
      axis.title = ggplot2::element_text(size=16)
    )

    graphics::hist(log10(ave_counts), breaks=100, main="", col="grey80",
                   xlab=expression(Log[10]~"average count"))
    graphics::abline(v=log10(min_ave_expression), col="blue", lwd=2, lty=2)
    normalisation_plots$initial_gene_filter <- grDevices::recordPlot()

    print(scater::plotQC(sce, type = "highest-expression", n=50) + fontsize)
    normalisation_plots$top_genes_qc <- grDevices::recordPlot()
  }

  numcells <- scater::nexprs(sce, byrow=TRUE)
  alt_keep <- numcells >= ncol(sce) * expressed_in_n_cells

  if (verbose) {
    graphics::smoothScatter(log10(ave_counts), numcells, xlab=expression(Log[10]~"average count"), ylab="Number of expressing cells")
    if (has_spike) {
      is_ercc <- SingleCellExperiment::isSpike(sce, type="ERCC")
      graphics::points(log10(ave_counts[is_ercc]), numcells[is_ercc], col="red", pch=16, cex=0.5)
    }

    normalisation_plots$cell_filtering <- grDevices::recordPlot()
  }

  sce <- sce[keep,]

  if (verbose) {
    normalisation_steps <- normalisation_steps %>%
      add_row(type = "gene_expression_filtering", ngenes = nrow(sce), ncells = ncol(sce))
    print(glue::glue("Gene filter: Genes - {nrow(sce)} Cells - {ncol(sce)}"))
  }

  ########################################
  # Normalise
  ########################################

  if (ncol(sce) >= 100) {
    # sizes <- c(20, 40, 60, 80)
    sizes <- ncol(sce) / 10
  } else {
    sizes <- ncol(sce)
  }

  sce <- scran::computeSumFactors(sce, sizes = sizes, positive = TRUE)
  sce <- sce[, sizeFactors(sce) > 0] # as mentioned in the scran documentation, ensure that size factors are higher than 0

  if (verbose) {
    graphics::plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
                   ylab="Library size (millions)", xlab="Size factor")
    normalisation_plots$size_factor <- grDevices::recordPlot()
  }

  if(has_spike) {
    sce <- scran::computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
    if(any(is.na(sizeFactors(sce, type="ERCC")))) {
      warning("Some cells do not have any spike-ins, this will cause an error further away. Remove spike-ins.")
    }
  }

  sce <- scater::normalise(sce)

  if (verbose) {
    normalisation_steps <- normalisation_steps %>%
      add_row(type = "normalisation", ngenes = nrow(sce), ncells = ncol(sce))
    print(glue::glue("Normalised: Genes - {nrow(sce)} Cells - {ncol(sce)}"))
  }

  ########################################
  # Select highly variable genes
  ########################################

  if (filter_hvg) {
    var_fit <- scran::trendVar(sce, method="spline", use.spikes=has_spike) # requires aroma.light
    var_out <- scran::decomposeVar(sce, var_fit)

    if (verbose) {
      graphics::plot(var_out$mean, var_out$total, pch=16, cex=0.6, xlab="Mean log-expression",
                     ylab="Variance of log-expression")
      o <- order(var_out$mean)
      graphics::lines(var_out$mean[o], var_out$tech[o], col="dodgerblue", lwd=2)

      if (has_spike) {
        cur_spike <- SingleCellExperiment::isSpike(sce)
        graphics::points(var_out$mean[cur_spike], var_out$total[cur_spike], col="red", pch=16)
      }

      normalisation_plots$gene_variance <- grDevices::recordPlot()

      normalisation_plots$gene_selection <- var_out %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes_string("FDR", "bio")) +
        ggplot2::geom_hline(yintercept = hvg_bio) +
        ggplot2::geom_vline(xintercept = hvg_fdr)
    }

    var_out <- var_out[order(var_out$bio, decreasing=TRUE),]
    hvg_out <- var_out[which(var_out$FDR <= hvg_fdr & var_out$bio >= hvg_bio),]
    if(nrow(hvg_out) < min_variable_fraction * nrow(var_out)) {
      hvg_out <- var_out[seq(1, ceiling(min_variable_fraction * nrow(var_out))), ]
    }

    # errors
    # if (verbose & nrow(hvg_out) >= 10) {
    #   normalisation_plots$top_genes <- scater::plotExpression(sce, rownames(hvg_out)[1:10]) + fontsize
    #   normalisation_plots$bottom_genes <- scater::plotExpression(sce, rownames(hvg_out)[(nrow(hvg_out)-10):nrow(hvg_out)]) + fontsize
    # }
    sce <- sce[rownames(hvg_out),]

    if (verbose) {
      normalisation_steps <- normalisation_steps %>%
        add_row(type = "gene_variability_filtering", ngenes = nrow(sce), ncells = ncol(sce))
      print(glue::glue("Variable genes filtered: Genes - {nrow(sce)} Cells - {ncol(sce)}"))
    }
  }

  expr_norm_filt <- Biobase::exprs(sce) %>% Matrix::t()
  counts_filt <- counts[rownames(expr_norm_filt),colnames(expr_norm_filt)]

  ########################################
  # Iterative filtering on variability
  ########################################
  repeat {
    gene_sds <- counts_filt %>% apply(2, stats::sd)
    cell_sds <- counts_filt %>% apply(1, stats::sd)

    genes_filtered <- which(gene_sds > 0, useNames = TRUE)
    cells_filtered <- which(cell_sds > 0, useNames = TRUE)
    expr_norm_filt <- expr_norm_filt[cells_filtered, genes_filtered]
    counts_filt <- counts_filt[cells_filtered, genes_filtered]

    if((min(c(gene_sds, 1), na.rm = TRUE) > 0) && (min(c(cell_sds, 1), na.rm = TRUE) > 0)){
      break
    }
  }

  if (verbose) {
    normalisation_steps <- normalisation_steps %>%
      add_row(type = "final_filtering", ngenes = ncol(expr_norm_filt), ncells = nrow(expr_norm_filt))
    print(glue::glue("Final filtering: Genes - {ncol(expr_norm_filt)} Cells - {nrow(expr_norm_filt)}"))
  }

  ########################################
  # Output
  ########################################

  if(verbose) {
    type <- NULL # satisfy r cmd check

    normalisation_plots$n_retained <-
      normalisation_steps %>%
      mutate(type = factor(type, levels=rev(type))) %>%
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
