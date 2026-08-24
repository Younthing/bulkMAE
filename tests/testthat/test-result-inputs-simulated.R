test_that("design and contrast helpers use aligned MAE metadata", {
  skip_if_not_installed("limma")
  mae <- make_simulated_mae(seed = 812L)

  design <- de_design(mae, "rna", ~ batch + condition)
  contrast <- de_contrast(
    mae,
    "rna",
    ~ batch + condition,
    contrasts = "conditiontreated"
  )

  expect_identical(rownames(design), rownames(mae_samples(mae, "rna")))
  expect_identical(rownames(contrast), colnames(design))
  expect_identical(colnames(contrast), "conditiontreated")
  expect_error(
    de_contrast(mae, "rna", ~ condition, character()),
    "contrast expressions"
  )
})

test_that("maSigPro design is generated only from aligned sample metadata", {
  mae <- make_simulated_mae(seed = 813L)
  design <- de_masigpro_design(
    mae,
    "rna",
    time = "time",
    replicate = "subject_id",
    group = "condition"
  )

  expect_identical(rownames(design), rownames(mae_samples(mae, "rna")))
  expect_identical(names(design), c("Time", "Replicates", "control", "treated"))
  expect_true(all(rowSums(design[, c("control", "treated")]) == 1L))
  expect_identical(
    unname(attr(design, "group_levels")),
    c("control", "treated")
  )

  binary <- mae_add_sample_data(
    mae,
    "rna",
    data.frame(
      control = as.integer(mae_samples(mae, "rna")$condition == "control"),
      treated = as.integer(mae_samples(mae, "rna")$condition == "treated"),
      row.names = rownames(mae_samples(mae, "rna"))
    )
  )
  binary_design <- de_masigpro_design(
    binary,
    "rna",
    time = "time",
    replicate = "subject_id",
    group = c("control", "treated")
  )
  expect_identical(
    as.matrix(binary_design[, c("control", "treated")]),
    as.matrix(design[, c("control", "treated")])
  )

  samples <- mae_samples(mae, "rna")
  samples$condition <- factor(
    rep(c("Time", "other"), length.out = nrow(samples))
  )
  renamed <- mae_from_matrix(
    mae_pull_assay(mae, "rna", "counts"),
    samples = samples,
    experiment = "rna"
  )
  reserved_design <- de_masigpro_design(
    renamed,
    "rna",
    time = "time",
    replicate = "subject_id",
    group = "condition"
  )
  expect_identical(anyDuplicated(names(reserved_design)), 0L)
  expect_identical(names(reserved_design)[1:3], c("Time", "Replicates", "Time.1"))
})

test_that("generic DE tables produce ranks and a complete selection vector", {
  native <- data.frame(
    log2FoldChange = c(2, -1.5, 0.1, NA_real_),
    lfcSE = c(0.5, 0.5, 0.5, NA_real_),
    stat = c(4, -3, 0.2, NA_real_),
    pvalue = c(1e-4, 0.002, 0.8, NA_real_),
    padj = c(4e-4, 0.004, 0.8, NA_real_),
    row.names = paste0("gene", 1:4)
  )

  table <- de_table(native)
  ranks <- de_ranks(native)
  selected <- de_selected(native, fdr = 0.01, min_abs_effect = 1)

  expect_identical(
    names(table),
    c(
      "feature_id", "effect", "standard_error", "statistic", "p_value",
      "adjusted_p_value"
    )
  )
  expect_identical(names(ranks), c("gene1", "gene3", "gene2"))
  expect_identical(
    selected,
    stats::setNames(c(TRUE, TRUE, FALSE, FALSE), paste0("gene", 1:4))
  )
  expect_error(de_selected(native[, 1:3]), "Adjusted p-values")
  invalid_se <- native
  invalid_se$lfcSE[[1L]] <- -0.5
  expect_error(de_table(invalid_se), "standard_error.*strictly positive")
  expect_error(
    de_table(data.frame(logFC = 1, F = -1, PValue = 0.5, row.names = "g1")),
    "F/LR statistics cannot be negative"
  )
  generic_f <- data.frame(
    logFC = c(1, -1),
    F = c(9, 4),
    PValue = c(0.01, 0.05),
    row.names = c("g1", "g2")
  )
  generic_f_table <- de_table(generic_f)
  expect_identical(generic_f_table$statistic, c(9, 4))
  expect_identical(
    attr(generic_f_table, "statistic_type"),
    "nondirectional_f_or_lr"
  )
  expect_error(de_ranks(generic_f), "Generic F/LR statistics.*non-directional")
  expect_identical(names(de_ranks(generic_f, column = "effect")), c("g1", "g2"))
})

test_that("limma results flow through de_table and enrichment ranks", {
  skip_if_not_installed("limma")
  mae <- make_simulated_mae(seed = 814L)
  fit <- de_limma(
    mae,
    "rna",
    formula = ~ condition,
    assay = "log_expression"
  )
  table <- de_table(fit, coef = "conditiontreated")
  ranks <- de_ranks(fit, coef = "conditiontreated")

  expect_identical(nrow(table), 1200L)
  expect_identical(names(ranks), table$feature_id[order(table$statistic, decreasing = TRUE)])
  expect_true(all(is.finite(table$standard_error)))
  expect_true(all(table$standard_error > 0))
  expect_error(de_table(fit), "coef.*required")
})

test_that("DESeq2 LRT statistics cannot be silently used as directional ranks", {
  skip_if_not_installed("DESeq2")
  mae <- mae_simulate(n_features = 80L, n_samples = 12L, seed = 816L)
  fit <- suppressMessages(de_deseq2(
    mae,
    "rna",
    design = ~ condition,
    test = "LRT",
    reduced = ~ 1
  ))
  result <- suppressMessages(de_deseq2_results(fit))
  table <- de_table(result)

  expect_identical(
    attr(table, "statistic_type"),
    "nondirectional_lrt"
  )
  expect_true(all(table$statistic[is.finite(table$statistic)] >= 0))
  expect_error(de_ranks(result), "LRT statistics are non-directional")
  effect_ranks <- de_ranks(result, column = "effect")
  expect_true(any(effect_ranks > 0))
  expect_true(any(effect_ranks < 0))
})

test_that("edgeR test statistics retain direction for downstream ranking", {
  skip_if_not_installed("edgeR")
  mae <- make_simulated_mae(seed = 815L)
  test <- suppressWarnings(de_edger(
    mae,
    "rna",
    formula = ~ condition,
    coef = "conditiontreated"
  ))
  table <- de_table(test)
  ranks <- de_ranks(test)

  expect_identical(nrow(table), nrow(test$table))
  expect_true(all(sign(table$statistic[table$effect != 0]) == sign(table$effect[table$effect != 0])))
  expect_identical(sort(names(ranks)), sort(table$feature_id[is.finite(table$statistic)]))
  expect_true(all(is.na(table$standard_error)))
  expect_error(
    meta_collect(
      list(cohort = test),
      feature = table$feature_id[[1L]]
    ),
    "lacks a finite positive standard error"
  )
  expect_true(any(de_selected(test, fdr = 0.2)))

  joint <- suppressWarnings(de_edger(
    mae,
    "rna",
    formula = ~ batch + condition,
    coef = c("batchbatch_b", "conditiontreated")
  ))
  expect_error(
    de_table(joint),
    "Multi-coefficient edgeR tests do not have one effect direction"
  )
  top_joint <- edgeR::topTags(joint, n = Inf)
  expect_error(
    de_table(top_joint, effect_column = "logFC.batchbatch_b"),
    "Multi-coefficient edgeR tests do not have one effect direction"
  )
})
