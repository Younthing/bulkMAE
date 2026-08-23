test_that("filter_expr removes simulated unexpressed features", {
  skip_if_not_installed("edgeR")

  source <- make_simulated_mae()
  counts <- mae_pull_assay(source, "rna", "counts")
  counts[seq_len(25L), ] <- 0L
  samples <- mae_samples(source, "rna")
  experiment <- mae_create_experiment(
    counts,
    data.frame(row.names = colnames(counts)),
    assay_name = "counts"
  )
  mae <- mae_create(list(rna = experiment), samples)

  keep <- filter_expr(
    mae,
    "rna",
    group = "condition",
    min_count = 10,
    min_total_count = 15
  )
  expect_type(keep, "logical")
  expect_length(keep, nrow(counts))
  expect_identical(names(keep), rownames(counts))
  expect_false(any(keep[seq_len(25L)]))
  expect_true(any(keep[-seq_len(25L)]))

  design <- stats::model.matrix(~ condition, data = samples)
  keep_from_design <- filter_expr(mae, "rna", design = design)
  expect_identical(names(keep_from_design), rownames(counts))
  expect_false(any(keep_from_design[seq_len(25L)]))
})

test_that("filter_expr rejects ambiguous grouping and invalid count assays", {
  skip_if_not_installed("edgeR")

  mae <- make_simulated_mae()
  samples <- mae_samples(mae, "rna")
  design <- stats::model.matrix(~ condition, data = samples)
  expect_error(
    filter_expr(mae, "rna", group = "condition", design = design),
    "only one"
  )
  expect_error(
    filter_expr(mae, "rna", group = rep("one", 2L)),
    "one value per sample"
  )

  counts <- mae_pull_assay(mae, "rna", "counts")
  counts[1L, 1L] <- NA_integer_
  invalid_experiment <- mae_create_experiment(
    counts,
    data.frame(row.names = colnames(counts)),
    assay_name = "counts"
  )
  invalid <- mae_create(list(rna = invalid_experiment), samples)
  expect_error(filter_expr(invalid, "rna"), "finite and non-negative")
})

test_that("qc_library reproduces direct simulated library summaries", {
  mae <- make_simulated_mae()
  counts <- mae_pull_assay(mae, "rna", "counts")
  observed <- qc_library(mae, "rna", "counts", detected_above = 5)

  expect_identical(
    names(observed),
    c("sample", "library_size", "detected_features", "zero_fraction")
  )
  expect_identical(observed$sample, colnames(counts))
  expect_equal(observed$library_size, unname(colSums(counts)))
  expect_equal(observed$detected_features, unname(colSums(counts > 5)))
  expect_equal(observed$zero_fraction, unname(colMeans(counts == 0)))
})

test_that("qc_library rejects non-finite assays and thresholds", {
  mae <- make_simulated_mae()
  samples <- mae_samples(mae, "rna")
  values <- mae_pull_assay(mae, "rna", "log_expression")
  values[1L, 1L] <- Inf
  experiment <- mae_create_experiment(
    values,
    data.frame(row.names = colnames(values)),
    assay_name = "values"
  )
  invalid <- mae_create(list(rna = experiment), samples)

  expect_error(qc_library(invalid, "rna", "values"), "non-finite")
  expect_error(qc_library(mae, "rna", detected_above = NA_real_), "finite number")
  expect_error(qc_library(mae, "rna", detected_above = c(0, 1)), "finite number")
})

test_that("qc_correlation returns an aligned sample correlation matrix", {
  mae <- make_simulated_mae()
  features <- sprintf("gene%04d", 1:150)
  expression <- mae_pull_assay(mae, "rna", "log_expression")[
    features,
    ,
    drop = FALSE
  ]

  observed <- qc_correlation(
    mae,
    "rna",
    "log_expression",
    method = "spearman",
    use = "everything",
    features = features
  )
  expect_identical(dim(observed), c(24L, 24L))
  expect_identical(rownames(observed), colnames(expression))
  expect_identical(colnames(observed), colnames(expression))
  expect_equal(observed, stats::cor(expression, method = "spearman"))
  expect_equal(unname(diag(observed)), rep(1, 24L))
  expect_equal(observed, t(observed))
})

test_that("qc_correlation validates feature selection and correlation inputs", {
  mae <- make_simulated_mae()

  expect_error(
    qc_correlation(
      mae,
      "rna",
      "log_expression",
      features = "gene0001"
    ),
    "At least two features"
  )
  expect_error(
    qc_correlation(
      mae,
      "rna",
      "log_expression",
      features = c("gene0001", "missing_gene")
    ),
    "Unknown features"
  )
  expect_error(
    qc_correlation(mae, "rna", "log_expression", method = "invalid"),
    "arg|method"
  )
})

test_that("reduce_pca embeds samples using selected variable features", {
  mae <- make_simulated_mae()
  result <- reduce_pca(
    mae,
    "rna",
    "log_expression",
    top_n = 75L,
    center = TRUE,
    scale. = TRUE
  )

  expect_s3_class(result, "prcomp")
  expect_identical(rownames(result$x), rownames(mae_samples(mae, "rna")))
  expect_identical(nrow(result$rotation), 75L)
  expect_true(all(is.finite(result$x)))
})

test_that("reduce_pca rejects invalid selection and non-finite values", {
  mae <- make_simulated_mae()
  expect_error(
    reduce_pca(mae, "rna", "log_expression", top_n = 0L),
    "positive integer"
  )

  samples <- mae_samples(mae, "rna")
  values <- mae_pull_assay(mae, "rna", "log_expression")
  values[1L, 1L] <- NA_real_
  experiment <- mae_create_experiment(
    values,
    data.frame(row.names = colnames(values)),
    assay_name = "values"
  )
  invalid <- mae_create(list(rna = experiment), samples)
  expect_error(reduce_pca(invalid, "rna", "values"), "finite assay")
})

test_that("reduce_mds returns a sample-aligned limma MDS object", {
  skip_if_not_installed("limma")

  mae <- make_simulated_mae()
  samples <- rownames(mae_samples(mae, "rna"))
  result <- reduce_mds(
    mae,
    "rna",
    "log_expression",
    features = sprintf("gene%04d", 1:200),
    top = 100L
  )

  expect_s4_class(result, "MDS")
  expect_identical(dim(result$distance.matrix.squared), c(24L, 24L))
  expect_identical(rownames(result$distance.matrix.squared), samples)
  expect_identical(colnames(result$distance.matrix.squared), samples)
  expect_true(all(is.finite(result$x)))
  expect_true(all(is.finite(result$y)))

  expect_error(
    reduce_mds(
      mae,
      "rna",
      "log_expression",
      features = "gene0001"
    ),
    "at least two features"
  )
})

test_that("reduce_umap is reproducible and labels simulated samples", {
  skip_if_not_installed("uwot")

  mae <- make_simulated_mae()
  arguments <- list(
    x = mae,
    experiment = "rna",
    assay = "log_expression",
    top_n = 80L,
    neighbors = 5L,
    components = 2L,
    seed = 17L,
    n_threads = 1L,
    n_sgd_threads = 1L,
    verbose = FALSE
  )
  first <- do.call(reduce_umap, arguments)
  second <- do.call(reduce_umap, arguments)

  expect_true(is.matrix(first))
  expect_identical(dim(first), c(24L, 2L))
  expect_identical(rownames(first), rownames(mae_samples(mae, "rna")))
  expect_equal(first, second, tolerance = 1e-10)
  expect_true(all(is.finite(first)))

  expect_error(
    reduce_umap(mae, "rna", "log_expression", neighbors = 24L),
    "samples minus one"
  )
  expect_error(
    reduce_umap(mae, "rna", "log_expression", components = 24L),
    "smaller than the sample count"
  )
})

test_that("reduce_tsne is reproducible and labels simulated samples", {
  skip_if_not_installed("Rtsne")

  mae <- make_simulated_mae()
  arguments <- list(
    x = mae,
    experiment = "rna",
    assay = "log_expression",
    top_n = 80L,
    dimensions = 2L,
    perplexity = 5,
    seed = 23L,
    max_iter = 250L,
    verbose = FALSE
  )
  first <- do.call(reduce_tsne, arguments)
  second <- do.call(reduce_tsne, arguments)

  expect_s3_class(first, "Rtsne")
  expect_identical(dim(first$Y), c(24L, 2L))
  expect_identical(rownames(first$Y), rownames(mae_samples(mae, "rna")))
  expect_equal(first$Y, second$Y, tolerance = 1e-10)
  expect_true(all(is.finite(first$Y)))

  expect_error(
    reduce_tsne(mae, "rna", "log_expression", perplexity = 8),
    "3 \\* perplexity"
  )
  expect_error(
    reduce_tsne(mae, "rna", "log_expression", dimensions = 24L),
    "smaller than the sample count"
  )
})

test_that("qc_outliers summarizes robust distances from simulated PCA", {
  mae <- make_simulated_mae()
  pca <- reduce_pca(mae, "rna", "log_expression", top_n = 100L)
  result <- qc_outliers(pca, components = 1:3, probability = 0.975)

  expect_identical(
    names(result),
    c("sample", "robust_distance", "cutoff", "flagged")
  )
  expect_identical(result$sample, rownames(pca$x))
  expect_length(result$robust_distance, 24L)
  expect_true(all(is.finite(result$robust_distance)))
  expect_equal(
    unique(result$cutoff),
    stats::qchisq(0.975, df = 3L)
  )
  expect_type(result$flagged, "logical")
})

test_that("qc_outliers rejects malformed PCA inputs and cutoffs", {
  mae <- make_simulated_mae()
  pca <- reduce_pca(mae, "rna", "log_expression", top_n = 50L)

  expect_error(qc_outliers(list()), "prcomp")
  expect_error(qc_outliers(pca, components = 999L), "unavailable")
  expect_error(qc_outliers(pca, components = integer()), "unavailable")
  expect_error(qc_outliers(pca, probability = 0), "between zero and one")
  expect_error(qc_outliers(pca, probability = 1), "between zero and one")
})
