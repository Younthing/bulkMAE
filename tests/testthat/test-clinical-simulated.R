test_that("the shared simulated MAE is reproducible without changing RNG state", {
  set.seed(481)
  caller_seed <- .Random.seed

  first <- make_simulated_mae(seed = 17L)
  expect_identical(.Random.seed, caller_seed)
  second <- make_simulated_mae(seed = 17L)

  expect_identical(.Random.seed, caller_seed)
  expect_identical(
    mae_pull_assay(first, "rna", "counts"),
    mae_pull_assay(second, "rna", "counts")
  )
  expect_identical(dim(mae_pull_assay(first, "rna", "counts")), c(1200L, 24L))
  expect_identical(
    mae_assays(first, "rna"),
    c("counts", "log_expression", "tpm", "weights")
  )
  expect_true(all(c(
    "condition", "batch", "subject_id", "time", "event",
    "risk_score", "outcome", "age"
  ) %in% names(mae_samples(first, "rna"))))
})

test_that("the shared fixture does not create a caller RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  make_simulated_mae(seed = 19L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("score_signature matches an explicit weighted sum", {
  mae <- make_simulated_mae()
  weights <- c(gene0001 = 1, gene0002 = -0.5, gene0003 = 0.25)
  expression <- mae_pull_assay(mae, "rna", "log_expression")

  score <- score_signature(
    mae,
    "rna",
    weights = weights,
    assay = "log_expression"
  )
  expected <- colSums(expression[names(weights), , drop = FALSE] * weights)

  expect_identical(names(score), colnames(expression))
  expect_equal(unname(score), unname(expected))
  expect_error(
    score_signature(
      mae,
      "rna",
      weights = c(gene0001 = 1, missing_gene = 1),
      assay = "log_expression"
    ),
    "Unknown features"
  )
  expect_error(
    score_signature(mae, "rna", weights = c(1, 2), assay = "log_expression"),
    "named numeric vector"
  )
})

test_that("survival wrappers fit native Kaplan-Meier and Cox models", {
  skip_if_not_installed("survival")
  mae <- make_simulated_mae()

  km <- surv_km(
    mae,
    "rna",
    survival::Surv(time, event) ~ condition
  )
  cox <- surv_cox(
    mae,
    "rna",
    survival::Surv(time, event) ~ age + gene0001,
    assay = "log_expression",
    features = "gene0001"
  )

  expect_s3_class(km, "survfit")
  expect_s3_class(cox, "coxph")
  expect_identical(sum(unname(km$n)), 24L)
  expect_true(all(is.finite(stats::coef(cox))))
  expect_error(surv_km(mae, "rna", "time ~ condition"), "formula")
  expect_error(
    surv_cox(
      mae,
      "rna",
      survival::Surv(time, event) ~ gene0001,
      features = "gene0001"
    ),
    "assay.*required"
  )
})

test_that("surv_roc returns a native timeROC result", {
  skip_if_not_installed("timeROC")
  mae <- make_simulated_mae()
  data <- mae_samples(mae, "rna")
  evaluation_times <- as.numeric(stats::quantile(data$time, c(0.25, 0.5)))

  result <- surv_roc(
    mae,
    "rna",
    time = "time",
    event = "event",
    marker = "risk_score",
    times = evaluation_times,
    cause = 1
  )

  expect_s3_class(result, "ipcwsurvivalROC")
  expect_true(all(is.finite(result$AUC)))
  expect_error(
    surv_roc(
      mae,
      "rna",
      time = "event",
      event = "event",
      marker = "risk_score",
      times = evaluation_times
    ),
    "must name different"
  )
})

test_that("penalized survival modelling is reproducible and validates folds", {
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  mae <- make_simulated_mae()
  features <- sprintf("gene%04d", seq_len(40L))

  set.seed(732)
  caller_seed <- .Random.seed
  first <- surv_penalized(
    mae,
    "rna",
    time = "time",
    event = "event",
    assay = "log_expression",
    features = features,
    folds = 3L,
    seed = 81L
  )
  second <- surv_penalized(
    mae,
    "rna",
    time = "time",
    event = "event",
    assay = "log_expression",
    features = features,
    folds = 3L,
    seed = 81L
  )

  expect_s3_class(first, "cv.glmnet")
  expect_identical(.Random.seed, caller_seed)
  expect_equal(first$lambda, second$lambda)
  expect_equal(first$cvm, second$cvm)
  expect_error(
    surv_penalized(
      mae,
      "rna",
      time = "time",
      event = "event",
      assay = "log_expression",
      features = features,
      folds = 17L
    ),
    "folds"
  )
})

test_that("ml_glmnet fits a name-aligned binary classifier", {
  skip_if_not_installed("glmnet")
  mae <- make_simulated_mae()
  features <- sprintf("gene%04d", seq_len(40L))

  fit <- ml_glmnet(
    mae,
    "rna",
    outcome = "outcome",
    assay = "log_expression",
    features = features,
    family = "binomial",
    folds = 3L,
    seed = 91L,
    foldid = rep(rep(seq_len(3L), each = 4L), times = 2L)
  )

  expect_s3_class(fit, "cv.glmnet")
  expect_true(all(is.finite(fit$cvm)))
  expect_error(
    ml_glmnet(
      mae,
      "rna",
      outcome = "outcome",
      assay = "log_expression",
      features = features,
      family = "binomial",
      folds = 13L
    ),
    "folds"
  )
})

test_that("meta_effect fits native models and rejects invalid standard errors", {
  skip_if_not_installed("metafor")
  effects <- c(0.25, 0.31, 0.12, 0.44, 0.18, 0.39, 0.21, 0.35)
  standard_errors <- c(0.10, 0.12, 0.11, 0.15, 0.09, 0.14, 0.10, 0.13)

  fit <- meta_effect(effects, standard_errors, method = "REML")

  expect_s3_class(fit, "rma.uni")
  expect_true(is.finite(as.numeric(fit$b)))
  expect_error(meta_effect(effects, standard_errors[-1L]), "equal lengths")
  expect_error(meta_effect(effects, replace(standard_errors, 1L, 0)), "strictly positive")
})

test_that("drug_query selects the strongest disjoint up and down genes", {
  statistics <- c(
    gene1 = 4,
    gene2 = 1,
    gene3 = -0.5,
    gene4 = -3,
    gene5 = 2,
    gene6 = -2
  )

  query <- drug_query(statistics, n = 2L)

  expect_identical(query$upset, c("gene1", "gene5"))
  expect_identical(query$downset, c("gene4", "gene6"))
  expect_length(intersect(query$upset, query$downset), 0L)
  expect_error(drug_query(c(gene1 = 2, gene2 = 1)), "positive and one negative")
})
