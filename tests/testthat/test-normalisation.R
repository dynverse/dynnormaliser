context("Normalisation")

tmpfile <- tempfile()
on.exit(unlink(tmpfile))

test_that("Testing normalise function", {
  num_genes <- 1000
  num_cells <- 1000

  gene_ids <- paste0("Gene", seq_len(num_genes))
  gene_ids[5:8] <- paste0("Mt-", seq_len(4))
  cell_ids <- paste0("Cell", seq_len(num_cells))

  counts <- matrix(round(2^rnorm(num_genes * num_cells, mean = 6, sd = 2)), ncol = num_genes, dimnames = list(cell_ids, gene_ids))
  counts[sample(c(T, F), length(counts), prob = c("T" = .1, "F" = .9), replace = TRUE)] <- 0
  counts[,1:4] <- counts[,1:4] + 10

  sink(tmpfile)
  pdf(tmpfile)
  normd <- normalise_filter_counts(
    counts = counts,
    verbose = TRUE,
    nmads = 3,
    min_ave_expression = 0.05,
    filter_hvg = TRUE,
    hvg_fdr = 0.05,
    hvg_bio = 0.5,
    min_variable_fraction = 0.1
  )
  dev.off()
  sink()

  expect_equal(dimnames(normd$counts), dimnames(normd$expression))
  expect_true(all(rownames(normd$counts) %in% rownames(counts)))
  expect_true(all(colnames(normd$counts) %in% colnames(counts)))
})



test_that("Also test for when there are no Mt genes", {
  num_genes <- 1000
  num_cells <- 1000

  gene_ids <- paste0("Gene", seq_len(num_genes))
  cell_ids <- paste0("Cell", seq_len(num_cells))

  counts <- matrix(round(2^rnorm(num_genes * num_cells, mean = 6, sd = 2)), ncol = num_genes, dimnames = list(cell_ids, gene_ids))
  counts[sample(c(T, F), length(counts), prob = c("T" = .1, "F" = .9), replace = TRUE)] <- 0
  counts[,1:4] <- counts[,1:4] + 10

  sink(tmpfile)
  pdf(tmpfile)
  normd <- normalise_filter_counts(
    counts = counts,
    verbose = TRUE,
    nmads = 3,
    min_ave_expression = 0.05,
    filter_hvg = TRUE,
    hvg_fdr = 0.05,
    hvg_bio = 0.5,
    min_variable_fraction = 0.1
  )
  dev.off()
  sink()

  expect_equal(dimnames(normd$counts), dimnames(normd$expression))
  expect_true(all(rownames(normd$counts) %in% rownames(counts)))
  expect_true(all(colnames(normd$counts) %in% colnames(counts)))
})
