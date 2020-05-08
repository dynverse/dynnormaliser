#' Normalisation
#'
#' @param counts The counts matrix, with features in columns
#' @export
normalise_filter_counts <- function(
  counts
) {
  as(log2(counts + 1), "dgCMatrix")
}
