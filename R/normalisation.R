#' Normalisation
#'
#' @param counts The counts matrix, with features in columns
#' @import Matrix
#' @export
normalise_filter_counts <- function(
  counts
) {
  list(
    counts = counts,
    expression = as(log2(counts + 1), "dgCMatrix")
  )
}
