.make_preprocessing_subset <- function(x, n_features = 200L) {
  sample_data <- mae_samples(x, "rna")
  feature_names <- rownames(mae_pull_assay(x, "rna", "counts"))[seq_len(n_features)]
  assay_names <- c("counts", "log_expression", "tpm", "weights")
  assays <- lapply(assay_names, function(assay) {
    mae_pull_assay(x, "rna", assay)[feature_names, , drop = FALSE]
  })
  names(assays) <- assay_names
  experiment <- mae_create_experiment(
    assays = assays,
    col_data = data.frame(row.names = rownames(sample_data))
  )
  mae_create(list(rna = experiment), sample_data)
}

test_that("edgeR calculates finite TMM normalization factors", {
  skip_if_not_installed("edgeR")
  mae <- .make_preprocessing_subset(make_simulated_mae())

  normalized <- normalize_tmm(mae, "rna", assay = "counts")

  expect_true(methods::is(normalized, "DGEList"))
  expect_identical(colnames(normalized), rownames(mae_samples(mae, "rna")))
  expect_true(all(is.finite(normalized$samples$norm.factors)))
  expect_true(all(normalized$samples$norm.factors > 0))
})

test_that("DESeq2 normalization and transformations preserve assay alignment", {
  skip_if_not_installed("DESeq2")
  mae <- make_simulated_mae()
  sample_names <- rownames(mae_samples(mae, "rna"))

  normalized <- normalize_deseq(mae, "rna", assay = "counts")
  vst <- transform_vst(
    mae,
    "rna",
    assay = "counts",
    blind = TRUE,
    fit_type = "mean"
  )
  rlog <- transform_rlog(
    mae,
    "rna",
    assay = "counts",
    blind = TRUE,
    fit_type = "mean"
  )

  expect_s4_class(normalized, "DESeqDataSet")
  expect_true(all(is.finite(DESeq2::sizeFactors(normalized))))
  expect_s4_class(vst, "DESeqTransform")
  expect_s4_class(rlog, "DESeqTransform")
  expect_identical(colnames(vst), sample_names)
  expect_identical(colnames(rlog), sample_names)
  expect_true(all(is.finite(SummarizedExperiment::assay(vst))))
  expect_true(all(is.finite(SummarizedExperiment::assay(rlog))))
})

test_that("voom returns aligned expression values and precision weights", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  mae <- .make_preprocessing_subset(make_simulated_mae())

  transformed <- transform_voom(
    mae,
    "rna",
    formula = ~ condition + batch,
    assay = "counts",
    plot = FALSE
  )

  expect_true(methods::is(transformed, "EList"))
  expect_identical(dim(transformed$E), c(200L, 24L))
  expect_identical(dim(transformed$weights), c(200L, 24L))
  expect_identical(colnames(transformed$E), rownames(mae_samples(mae, "rna")))
  expect_true(all(is.finite(transformed$E)))
  expect_true(all(is.finite(transformed$weights)))
})

test_that("continuous and count batch adjustment return aligned matrices", {
  suppressWarnings(skip_if_not_installed("sva"))
  skip_if_not_installed("limma")
  mae <- .make_preprocessing_subset(make_simulated_mae())
  expected_names <- colnames(mae_pull_assay(mae, "rna", "counts"))

  combat <- adjust_combat(
    mae,
    "rna",
    batch = "batch",
    assay = "log_expression",
    preserve = ~ condition,
    parametric = TRUE
  )
  removed <- adjust_batch(
    mae,
    "rna",
    batch = "batch",
    assay = "log_expression",
    preserve = ~ condition
  )
  combat_seq <- adjust_combatseq(
    mae,
    "rna",
    batch = "batch",
    assay = "counts",
    group = "condition"
  )

  for (result in list(combat, removed, combat_seq)) {
    expect_true(is.matrix(result))
    expect_identical(dim(result), c(200L, 24L))
    expect_identical(colnames(result), expected_names)
    expect_true(all(is.finite(result)))
  }
  expect_true(all(combat_seq >= 0))
})

test_that("SVA estimates a requested surrogate variable", {
  suppressWarnings(skip_if_not_installed("sva"))
  mae <- .make_preprocessing_subset(make_simulated_mae())

  result <- adjust_sva(
    mae,
    "rna",
    assay = "log_expression",
    full = ~ condition + batch,
    null = ~ batch,
    n_surrogates = 1L,
    method = "irw"
  )

  expect_true(is.list(result))
  expect_identical(result$n.sv, 1L)
  expect_length(result$sv, 24L)
  expect_true(all(is.finite(result$sv)))
})

test_that("RUVg accepts named negative-control features", {
  suppressWarnings(skip_if_not_installed("RUVSeq"))
  mae <- .make_preprocessing_subset(make_simulated_mae())
  controls <- sprintf("gene%04d", 151:180)

  result <- adjust_ruv(
    mae,
    "rna",
    controls = controls,
    k = 1L,
    assay = "counts"
  )

  expect_true(is.list(result))
  expect_identical(dim(result$W), c(24L, 1L))
  expect_identical(dim(result$normalizedCounts), c(200L, 24L))
  expect_identical(
    colnames(result$normalizedCounts),
    rownames(mae_samples(mae, "rna"))
  )
  expect_true(all(is.finite(result$W)))
})

test_that("variancePartition returns one aligned row per feature", {
  skip_if_not_installed("variancePartition")
  mae <- .make_preprocessing_subset(make_simulated_mae(), n_features = 60L)

  result <- de_variance(
    mae,
    "rna",
    formula = ~ age + risk_score,
    assay = "log_expression",
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(is.data.frame(result))
  expect_identical(rownames(result), sprintf("gene%04d", seq_len(60L)))
  expect_true(all(is.finite(as.matrix(result))))
})
