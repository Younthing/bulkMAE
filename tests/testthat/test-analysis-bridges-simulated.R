test_that("installed simulation entry point provides a complete aligned MAE", {
  set.seed(991)
  caller_seed <- .Random.seed
  mae <- mae_simulate(n_features = 80L, n_samples = 8L, seed = 91L)

  expect_s4_class(mae, "MultiAssayExperiment")
  expect_identical(mae_experiments(mae), "rna")
  expect_identical(mae_assays(mae, "rna"), c("counts", "log_expression", "tpm"))
  expect_identical(dim(mae_pull_assay(mae, "rna", "counts")), c(80L, 8L))
  expect_identical(rownames(mae_feature_data(mae, "rna")), sprintf("gene%04d", 1:80))
  expect_equal(
    colSums(mae_pull_assay(mae, "rna", "tpm")),
    stats::setNames(rep(1e6, 8L), sprintf("sample%02d", 1:8))
  )
  expect_identical(.Random.seed, caller_seed)
})

test_that("the default simulation contains a detectable time-course signal", {
  skip_if_not_installed("maSigPro")
  mae <- mae_simulate(seed = 1L)
  design <- de_masigpro_design(
    mae,
    "rna",
    time = "timepoint",
    replicate = "replicate",
    group = "condition"
  )
  invisible(utils::capture.output(
    result <- suppressMessages(suppressWarnings(de_masigpro(
      mae,
      "rna",
      edesign = design,
      assay = "log_expression"
    )))
  ))
  expect_type(result, "list")
  expect_gt(length(result), 0L)
})

test_that("survival formulas and derived sample scores stay inside the public API", {
  skip_if_not_installed("survival")
  mae <- mae_simulate(n_features = 80L, n_samples = 12L, seed = 92L)
  score <- score_signature(
    mae,
    "rna",
    weights = sprintf("gene%04d", 1:10),
    assay = "log_expression"
  )
  derived <- data.frame(
    signature_score = score,
    signature_group = factor(ifelse(score >= stats::median(score), "high", "low")),
    row.names = names(score)
  )
  augmented <- mae_add_sample_data(mae, "rna", derived)
  km_formula <- surv_formula(
    time = "survival_time",
    event = "event",
    predictors = "signature_group"
  )
  cox_formula <- surv_formula(
    time = "survival_time",
    event = "event",
    predictors = "signature_score"
  )

  km <- surv_km(augmented, "rna", km_formula)
  cox <- surv_cox(augmented, "rna", cox_formula)
  expect_s3_class(km, "survfit")
  expect_s3_class(cox, "coxph")
  expect_true("signature_score" %in% names(mae_samples(augmented, "rna")))
  expect_false("signature_score" %in% names(mae_samples(mae, "rna")))
})

test_that("multiple native DE fits flow into meta-analysis through public helpers", {
  skip_if_not_installed("limma")
  skip_if_not_installed("metafor")
  first <- mae_simulate(n_features = 80L, n_samples = 12L, seed = 93L)
  second <- mae_simulate(n_features = 80L, n_samples = 12L, seed = 94L)
  fits <- list(
    cohort_a = de_limma(first, "rna", ~ condition, assay = "log_expression"),
    cohort_b = de_limma(second, "rna", ~ condition, assay = "log_expression")
  )
  effects <- meta_collect(
    fits,
    feature = "gene0001",
    coef = "conditiontreated"
  )
  meta <- meta_effect(effects, method = "FE")

  expect_identical(effects$study, c("cohort_a", "cohort_b"))
  expect_true(all(is.finite(effects$effect)))
  expect_true(all(effects$standard_error > 0))
  expect_s3_class(meta, "rma")

  multivariable <- list(
    cohort_a = de_limma(first, "rna", ~ batch + condition, assay = "log_expression"),
    cohort_b = de_limma(second, "rna", ~ batch + condition, assay = "log_expression")
  )
  coefficients <- list(
    cohort_b = "batchbatch_b",
    cohort_a = "conditiontreated"
  )
  aligned <- meta_collect(
    multivariable,
    feature = "gene0001",
    coef = coefficients
  )
  expect_equal(
    aligned$effect,
    c(
      de_table(multivariable$cohort_a, coef = "conditiontreated")$effect[[1L]],
      de_table(multivariable$cohort_b, coef = "batchbatch_b")$effect[[1L]]
    )
  )
  expect_error(
    meta_collect(
      multivariable,
      feature = "gene0001",
      coef = list(wrong_a = 2L, wrong_b = 2L)
    ),
    "same unique names"
  )
})

test_that("activity, WGCNA, and NMF outputs have public downstream accessors", {
  activity <- expand.grid(
    source = c("pathway_a", "pathway_b"),
    condition = c("sample1", "sample2", "sample3"),
    statistic = c("ulm", "mlm"),
    stringsAsFactors = FALSE
  )
  activity$score <- seq_len(nrow(activity)) / 10
  matrix <- activity_matrix(activity, statistic = "ulm")
  expect_identical(dim(matrix), c(2L, 3L))
  expect_error(activity_matrix(activity), "statistic.*required")

  modules <- coexpr_modules(list(colors = c(gene1 = "blue", gene2 = "brown")))
  expect_identical(modules, c(gene1 = "blue", gene2 = "brown"))
  expect_error(coexpr_modules(list(colors = c("blue", "brown"))), "named")
})

test_that("NMF native fits expose named sample and feature classes", {
  skip_if_not_installed("NMF")
  mae <- mae_simulate(n_features = 30L, n_samples = 8L, seed = 95L)
  fit <- suppressMessages(cluster_nmf(
    mae,
    "rna",
    assay = "tpm",
    rank = 2L,
    runs = 1L,
    seed = 95L,
    .options = list(maxIter = 30L)
  ))

  sample_classes <- cluster_nmf_classes(fit, what = "samples")
  feature_classes <- cluster_nmf_classes(fit, what = "features")

  expect_identical(names(sample_classes), sprintf("sample%02d", 1:8))
  expect_identical(names(feature_classes), sprintf("gene%04d", 1:30))
  expect_true(all(as.integer(sample_classes) %in% 1:2))
  expect_true(all(as.integer(feature_classes) %in% 1:2))
})
