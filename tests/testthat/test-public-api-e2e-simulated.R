test_that("public API runs matrix to filtered DE and offline enrichment", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("fgsea")

  simulated <- mae_simulate(
    n_features = 180L,
    n_samples = 16L,
    seed = 1201L
  )
  mae <- mae_from_matrix(
    expression = mae_pull_assay(simulated, "rna", "counts"),
    samples = mae_samples(simulated, "rna"),
    row_data = mae_feature_data(simulated, "rna"),
    experiment = "rna",
    assay = "counts"
  )
  keep <- filter_expr(
    mae,
    "rna",
    assay = "counts",
    group = "condition"
  )
  filtered <- mae_subset_features(mae, "rna", names(keep)[keep])
  voom <- transform_voom(
    filtered,
    "rna",
    formula = ~ condition,
    assay = "counts",
    plot = FALSE
  )
  normalized <- mae_add_assay(
    filtered,
    "rna",
    value = voom,
    name = "voom"
  )
  design <- de_design(normalized, "rna", ~ condition)
  contrast <- suppressMessages(de_contrast(
    normalized,
    "rna",
    formula = ~ condition,
    contrasts = "conditiontreated"
  ))
  fit <- de_limma(
    normalized,
    "rna",
    formula = ~ condition,
    assay = "voom",
    contrasts = contrast
  )
  table <- de_table(fit)
  ranks <- de_ranks(fit)
  selected <- de_selected(fit, fdr = 1, min_abs_effect = 0)
  ranked_features <- names(ranks)
  gene_sets <- gene_sets_prepare(
    list(
      top = ranked_features[seq_len(30L)],
      bottom = tail(ranked_features, 30L),
      mixed = c(
        ranked_features[31:45],
        tail(ranked_features, 15L)
      )
    ),
    min_size = 15L,
    max_size = 40L
  )
  enrichment <- suppressWarnings(enrich_fgsea(
    ranks = ranks,
    pathways = gene_sets,
    min_size = 10L,
    max_size = 50L,
    eps = 0
  ))

  expect_true("voom" %in% mae_assays(normalized, "rna"))
  expect_identical(rownames(design), rownames(mae_samples(normalized, "rna")))
  expect_identical(rownames(contrast), colnames(design))
  expect_identical(nrow(table), sum(keep))
  expect_identical(sort(names(ranks)), sort(table$feature_id))
  expect_identical(names(selected), table$feature_id)
  expect_true(is.data.frame(enrichment))
  expect_setequal(enrichment$pathway, names(gene_sets))
})

test_that("public SVA covariates flow back into differential modelling", {
  suppressWarnings(skip_if_not_installed("sva"))
  skip_if_not_installed("limma")

  mae <- mae_simulate(
    n_features = 100L,
    n_samples = 12L,
    seed = 1202L
  )
  adjustment <- suppressMessages(adjust_sva(
    mae,
    "rna",
    assay = "log_expression",
    full = ~ condition + batch,
    null = ~ batch,
    n_surrogates = 1L
  ))
  covariates <- adjust_covariates(adjustment)
  augmented <- mae_add_sample_data(mae, "rna", covariates)
  fit <- de_limma(
    augmented,
    "rna",
    formula = ~ batch + SV1 + condition,
    assay = "log_expression"
  )
  table <- de_table(fit, coef = "conditiontreated")

  expect_identical(rownames(covariates), rownames(mae_samples(mae, "rna")))
  expect_false("SV1" %in% names(mae_samples(mae, "rna")))
  expect_true("SV1" %in% names(mae_samples(augmented, "rna")))
  expect_identical(nrow(table), 100L)
  expect_true(any(is.finite(table$statistic)))
})

test_that("public signature scoring flows into KM and Cox models", {
  skip_if_not_installed("survival")

  mae <- mae_simulate(
    n_features = 100L,
    n_samples = 20L,
    seed = 1203L
  )
  score <- score_signature(
    mae,
    "rna",
    weights = sprintf("gene%04d", seq_len(12L)),
    assay = "log_expression"
  )
  signature_data <- data.frame(
    signature_score = score,
    signature_group = factor(ifelse(
      score >= median(score),
      "high",
      "low"
    )),
    row.names = names(score),
    check.names = FALSE
  )
  augmented <- mae_add_sample_data(mae, "rna", signature_data)
  km_formula <- surv_formula(
    time = "survival_time",
    event = "event",
    predictors = "signature_group"
  )
  cox_formula <- surv_formula(
    time = "survival_time",
    event = "event",
    predictors = c("signature_score", "age")
  )
  km <- surv_km(augmented, "rna", km_formula)
  cox <- surv_cox(augmented, "rna", cox_formula)

  expect_identical(names(score), rownames(mae_samples(mae, "rna")))
  expect_true(all(c(
    "signature_score",
    "signature_group"
  ) %in% names(mae_samples(augmented, "rna"))))
  expect_s3_class(km, "survfit")
  expect_s3_class(cox, "coxph")
})

test_that("public metadata design runs a complete maSigPro workflow", {
  skip_if_not_installed("maSigPro")

  mae <- mae_simulate(
    n_features = 60L,
    n_samples = 18L,
    seed = 1204L
  )
  edesign <- de_masigpro_design(
    mae,
    "rna",
    time = "timepoint",
    replicate = "replicate",
    group = "condition"
  )
  result <- de_masigpro(
    mae,
    "rna",
    edesign = edesign,
    assay = "log_expression",
    degree = 2L,
    q = 1,
    counts = FALSE,
    alpha = 1,
    r_squared = 0,
    min_observations = 9L
  )

  expect_identical(rownames(edesign), rownames(mae_samples(mae, "rna")))
  expect_identical(
    names(edesign),
    c("Time", "Replicates", "control", "treated")
  )
  expect_type(result, "list")
  expect_gt(length(result), 0L)
})

test_that("public GSVA scores become an MAE experiment for limma", {
  skip_if_not_installed("GSVA", minimum_version = "2.0.0")
  skip_if_not_installed("limma")

  mae <- mae_simulate(
    n_features = 160L,
    n_samples = 12L,
    seed = 1205L
  )
  gene_sets <- gene_sets_prepare(
    list(
      signal = sprintf("gene%04d", seq_len(35L)),
      other = sprintf("gene%04d", 36:70),
      mixed = c(
        sprintf("gene%04d", seq_len(15L)),
        sprintf("gene%04d", 71:90)
      )
    ),
    min_size = 15L,
    max_size = 40L
  )
  scores <- score_gsva(
    mae,
    "rna",
    gene_sets = gene_sets,
    assay = "log_expression",
    kcdf = "Gaussian",
    min_size = 10L,
    max_size = 50L,
    verbose = FALSE
  )
  augmented <- mae_add_experiment(
    mae,
    value = scores,
    name = "pathway",
    assay_name = "score",
    source_experiment = "rna"
  )
  fit <- de_limma(
    augmented,
    "pathway",
    formula = ~ condition,
    assay = "score"
  )
  table <- de_table(fit, coef = "conditiontreated")

  expect_true("pathway" %in% mae_experiments(augmented))
  expect_identical(mae_assays(augmented, "pathway"), "score")
  expect_setequal(table$feature_id, names(gene_sets))
  expect_true(all(is.finite(table$statistic)))
})

test_that("public WGCNA modules flow into cohort preservation", {
  skip_if_not_installed("WGCNA")

  reference <- mae_simulate(
    n_features = 80L,
    n_samples = 16L,
    seed = 1301L
  )
  test_cohort <- mae_simulate(
    n_features = 80L,
    n_samples = 16L,
    seed = 1302L
  )
  fit <- suppressWarnings(suppressMessages(coexpr_wgcna(
    reference,
    "rna",
    assay = "log_expression",
    power = 2,
    min_module_size = 5L,
    merge_cut_height = 0.35,
    numeric_labels = FALSE,
    seed = 1303L,
    verbose = 0,
    maxBlockSize = 100L,
    pamRespectsDendro = FALSE
  )))
  modules <- coexpr_modules(fit)
  preservation <- suppressWarnings(suppressMessages(coexpr_preservation(
    reference,
    "rna",
    test_cohort,
    "rna",
    module_colors = modules,
    reference_assay = "log_expression",
    test_assay = "log_expression",
    permutations = 1L,
    network_type = "signed",
    seed = 1304L,
    verbose = 0,
    save_permuted_statistics = FALSE
  )))

  expect_named(modules)
  expect_setequal(
    names(modules),
    rownames(mae_pull_assay(reference, "rna", "log_expression"))
  )
  expect_gt(length(unique(modules)), 1L)
  expect_type(preservation, "list")
  expect_true("preservation" %in% names(preservation))
})
