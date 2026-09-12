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

