.make_differential_subset <- function(x, n_features = 300L) {
  sample_data <- mae_samples(x, "rna")
  sample_data$condition <- factor(sample_data$condition)
  sample_data$batch <- factor(sample_data$batch)
  sample_data$subject_id <- factor(sample_data$subject_id)
  feature_names <- rownames(mae_pull_assay(x, "rna", "counts"))[seq_len(n_features)]
  assays <- list(
    counts = mae_pull_assay(x, "rna", "counts")[feature_names, , drop = FALSE],
    log_expression = mae_pull_assay(
      x,
      "rna",
      "log_expression"
    )[feature_names, , drop = FALSE]
  )
  experiment <- mae_create_experiment(
    assays = assays,
    col_data = data.frame(row.names = rownames(sample_data))
  )
  mae_create(list(rna = experiment), sample_data)
}

.make_exon_fixture <- function(seed = 210L) {
  source <- .make_differential_subset(
    make_simulated_mae(seed = seed),
    n_features = 180L
  )
  sample_data <- mae_samples(source, "rna")
  exon_names <- sprintf("exon%03d", seq_len(180L))
  counts <- mae_pull_assay(source, "rna", "counts")
  expression <- mae_pull_assay(source, "rna", "log_expression")
  rownames(counts) <- exon_names
  rownames(expression) <- exon_names
  experiment <- mae_create_experiment(
    assays = list(counts = counts, log_expression = expression),
    col_data = data.frame(row.names = rownames(sample_data))
  )
  list(
    mae = mae_create(list(exon = experiment), sample_data),
    gene_id = stats::setNames(
      rep(sprintf("splice_gene%02d", seq_len(60L)), each = 3L),
      exon_names
    ),
    feature_id = stats::setNames(exon_names, exon_names)
  )
}

.make_timecourse_fixture <- function(seed = 310L) {
  source <- .make_differential_subset(
    make_simulated_mae(seed = seed),
    n_features = 120L
  )
  sample_names <- colnames(mae_pull_assay(source, "rna", "log_expression"))
  time <- rep(rep(0:3, each = 3L), times = 2L)
  replicates <- rep(seq_len(3L), times = 8L)
  group <- rep(c("control", "treated"), each = 12L)
  edesign <- data.frame(
    Time = time,
    Replicates = replicates,
    control = as.integer(group == "control"),
    treated = as.integer(group == "treated"),
    row.names = sample_names,
    check.names = FALSE
  )
  expression <- mae_pull_assay(source, "rna", "log_expression")
  expression[seq_len(30L), ] <- expression[seq_len(30L), ] +
    outer(rep(1.4, 30L), time * as.integer(group == "treated"))
  sample_data <- data.frame(
    condition = group,
    timepoint = time,
    replicate = replicates,
    row.names = sample_names,
    check.names = FALSE
  )
  experiment <- mae_create_experiment(
    assays = list(log_expression = expression),
    col_data = data.frame(row.names = sample_names)
  )
  list(
    mae = mae_create(list(rna = experiment), sample_data),
    edesign = edesign
  )
}

test_that("DESeq2 fitting and result extraction preserve feature alignment", {
  skip_if_not_installed("DESeq2")
  mae <- .make_differential_subset(make_simulated_mae())

  fit <- de_deseq2(
    mae,
    "rna",
    design = ~ batch + condition,
    assay = "counts",
    fitType = "mean",
    quiet = TRUE
  )
  result <- de_deseq2_results(
    fit,
    contrast = c("condition", "treated", "control")
  )

  expect_s4_class(fit, "DESeqDataSet")
  expect_s4_class(result, "DESeqResults")
  expect_identical(rownames(result), sprintf("gene%04d", seq_len(300L)))
  expect_true(any(is.finite(result$stat)))
})

test_that("edgeR runs both fit-only and coefficient-test paths", {
  skip_if_not_installed("edgeR")
  mae <- .make_differential_subset(make_simulated_mae())

  fit <- de_edger(
    mae,
    "rna",
    formula = ~ batch + condition,
    assay = "counts",
    filter = FALSE,
    robust = FALSE
  )
  tested <- de_edger(
    mae,
    "rna",
    formula = ~ batch + condition,
    assay = "counts",
    filter = FALSE,
    robust = FALSE,
    coef = "conditiontreated"
  )

  expect_true(methods::is(fit, "DGEGLM"))
  expect_true(methods::is(tested, "DGELRT"))
  expect_identical(rownames(fit$coefficients), sprintf("gene%04d", seq_len(300L)))
})

test_that("voom and continuous limma fits retain native result classes", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  mae <- .make_differential_subset(make_simulated_mae())

  voom_fit <- de_voom(
    mae,
    "rna",
    formula = ~ batch + condition,
    assay = "counts",
    filter = FALSE,
    voom_plot = FALSE
  )
  continuous_fit <- de_limma(
    mae,
    "rna",
    formula = ~ batch + condition,
    assay = "log_expression"
  )

  expect_true(methods::is(voom_fit, "MArrayLM"))
  expect_true(methods::is(continuous_fit, "MArrayLM"))
  expect_identical(rownames(voom_fit$coefficients), rownames(continuous_fit$coefficients))
  expect_true(all(is.finite(continuous_fit$coefficients)))
})

test_that("diffSplice runs on both limma and edgeR exon fits", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  fixture <- .make_exon_fixture()

  limma_fit <- de_limma(
    fixture$mae,
    "exon",
    formula = ~ batch + condition,
    assay = "log_expression"
  )
  edge_fit <- de_edger(
    fixture$mae,
    "exon",
    formula = ~ batch + condition,
    assay = "counts",
    filter = FALSE,
    robust = FALSE
  )
  limma_result <- dtu_diffsplice(
    limma_fit,
    gene_id = fixture$gene_id,
    feature_id = fixture$feature_id
  )
  edge_result <- dtu_diffsplice(
    edge_fit,
    gene_id = fixture$gene_id,
    feature_id = fixture$feature_id,
    coef = "conditiontreated"
  )

  expect_true(methods::is(limma_result, "MArrayLM"))
  expect_true(methods::is(edge_result, "DGELRT"))
  expect_identical(rownames(limma_result$coefficients), names(fixture$gene_id))
  expect_identical(as.character(edge_result$genes$ExonID), unname(fixture$feature_id))
  expect_identical(as.character(edge_result$genes$GeneID), unname(fixture$gene_id))
})

test_that("maSigPro runs a complete aligned time-course workflow", {
  skip_if_not_installed("maSigPro")
  fixture <- .make_timecourse_fixture()

  result <- de_masigpro(
    fixture$mae,
    "rna",
    edesign = fixture$edesign,
    assay = "log_expression",
    degree = 2L,
    q = 1,
    counts = FALSE,
    alpha = 1,
    r_squared = 0,
    min_observations = 12L
  )

  expect_true(is.list(result))
  expect_true(length(result) > 0L)
})

test_that("dream fits a repeated-measures model with native moderation", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("lme4")
  skip_if_not_installed("limma")
  skip_if_not_installed("variancePartition")
  mae <- .make_differential_subset(make_simulated_mae(), n_features = 300L)

  fit <- de_dream(
    mae,
    "rna",
    formula = ~ condition + (1 | subject_id),
    assay = "counts",
    filter = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  expect_true(methods::is(fit, "MArrayLM"))
  expect_identical(rownames(fit$coefficients), sprintf("gene%04d", seq_len(300L)))
  expect_true(all(is.finite(fit$coefficients)))
})
