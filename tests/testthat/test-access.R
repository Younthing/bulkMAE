test_that("the access contract preserves sample order", {
  mae <- make_toy_mae()
  expect_silent(mae_validate(mae, "rna"))
  expect_identical(mae_experiments(mae), "rna")
  expect_identical(
    mae_assays(mae, "rna"),
    c("counts", "log_expression")
  )
  expect_identical(
    colnames(mae_pull_assay(mae, "rna", "counts")),
    rownames(mae_samples(mae, "rna"))
  )
})

test_that("library metrics are deterministic", {
  metrics <- qc_library(make_toy_mae(), "rna")
  expect_named(
    metrics,
    c("sample", "library_size", "detected_features", "zero_fraction")
  )
  expect_equal(metrics$library_size, c(6, 15, 24, 33))
})

test_that("unknown experiments and assays fail clearly", {
  mae <- make_toy_mae()
  expect_error(mae_pull_experiment(mae, "missing"), "Unknown experiment")
  expect_error(mae_pull_assay(mae, "rna", "missing"), "Unknown assay")
})

test_that("explicit sampleMap aligns primary metadata", {
  counts <- matrix(
    seq_len(6),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), c("aliquot1", "aliquot2"))
  )
  aliquots <- data.frame(row.names = colnames(counts))
  primary <- data.frame(
    condition = c("control", "treated"),
    row.names = c("patient1", "patient2")
  )
  map <- data.frame(
    assay = "rna",
    primary = c("patient2", "patient1"),
    colname = c("aliquot1", "aliquot2")
  )
  se <- mae_create_experiment(counts, aliquots)
  mae <- mae_create(list(rna = se), primary, sample_map = map)

  expect_identical(
    mae_samples(mae, "rna")$condition,
    c("treated", "control")
  )
})

test_that("mae_variable_features ranks finite-variance rows", {
  mae <- make_toy_mae(n_features = 6L, n_samples = 4L)
  ranked <- matrix(
    c(
      1, 1, 1, 1,
      1, 2, 3, 4,
      10, 20, 30, 40,
      0, 0, 1, 1,
      5, 5, 5, 6,
      100, 0, 0, 0
    ),
    nrow = 6L,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:6), paste0("sample", 1:4))
  )
  mae <- mae_add_assay(mae, "rna", ranked, name = "ranked")
  features <- mae_variable_features(mae, "rna", "ranked", top_n = 3L)
  variance <- apply(ranked, 1L, stats::var)
  expect_identical(features, names(sort(variance[variance > 0], decreasing = TRUE))[1:3])
  expect_false("gene1" %in% features)
  expect_error(
    mae_variable_features(mae, "rna", "ranked", top_n = 0L),
    "positive integer"
  )
})

test_that("mae_deconv_toy is synthetic and sample-aligned", {
  samples <- paste0("sample", 1:4)
  genes <- paste0("ENSG", 1:12)
  toy <- mae_deconv_toy(samples, genes = genes, seed = 7L)

  expect_identical(rownames(toy$counts), genes)
  expect_identical(colnames(toy$fractions), samples)
  expect_identical(rownames(toy$cell_data), colnames(toy$counts))
  expect_true(all(c("cell_type", "sample_id") %in% names(toy$cell_data)))
  expect_equal(colSums(toy$fractions), rep(1, 4L), tolerance = 1e-8, ignore_attr = TRUE)
  expect_match(toy$source, "Synthetic")
  expect_error(mae_deconv_toy(c("a", "a")), "unique")
  expect_error(mae_deconv_toy(samples, n_donors = 1L), "at least 2")

  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    reference <- deconv_reference(
      toy$counts,
      toy$cell_data,
      cell_type = "cell_type",
      sample = "sample_id"
    )
    expect_true(methods::is(reference, "SingleCellExperiment"))
    expect_identical(dim(reference), c(12L, ncol(toy$counts)))
  }
})
