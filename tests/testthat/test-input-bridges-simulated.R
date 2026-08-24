.make_input_bridge_mae <- function() {
  features <- paste0("gene", seq_len(6L))
  samples <- paste0("sample", seq_len(4L))
  counts <- matrix(
    seq_len(24L),
    nrow = length(features),
    dimnames = list(features, samples)
  )
  storage.mode(counts) <- "integer"
  sample_data <- data.frame(
    condition = rep(c("control", "treated"), each = 2L),
    row.names = samples,
    check.names = FALSE
  )
  feature_data <- data.frame(
    gene_length = seq(600, 1100, length.out = length(features)),
    row.names = features,
    check.names = FALSE
  )
  mae_from_matrix(
    counts,
    sample_data,
    row_data = feature_data,
    experiment = "rna",
    assay = "counts"
  )
}

test_that("mae_from_matrix creates a complete aligned one-experiment MAE", {
  mae <- .make_input_bridge_mae()

  expect_s4_class(mae, "MultiAssayExperiment")
  expect_identical(mae_experiments(mae), "rna")
  expect_identical(mae_assays(mae, "rna"), "counts")
  expect_identical(
    colnames(mae_pull_assay(mae, "rna", "counts")),
    rownames(mae_samples(mae, "rna"))
  )
  expect_identical(
    mae_feature_data(mae, "rna")$gene_length,
    seq(600, 1100, length.out = 6L)
  )
  raw_leaf <- MultiAssayExperiment::experiments(mae)[["rna"]]
  expect_identical(names(SummarizedExperiment::colData(raw_leaf)), character())

  counts <- mae_pull_assay(mae, "rna", "counts")
  reversed <- mae_samples(mae, "rna")[rev(colnames(counts)), , drop = FALSE]
  expect_error(
    mae_from_matrix(counts, reversed),
    "exactly match"
  )
})

test_that("feature metadata bridges align names and protect existing columns", {
  mae <- .make_input_bridge_mae()
  features <- rownames(mae_pull_assay(mae, "rna", "counts"))
  symbols <- stats::setNames(toupper(rev(features)), rev(features))

  with_symbols <- mae_add_feature_data(
    mae,
    "rna",
    value = symbols,
    name = "symbol"
  )
  observed <- mae_feature_data(with_symbols, "rna")
  expect_identical(rownames(observed), features)
  expect_identical(observed$symbol, toupper(features))
  expect_false("symbol" %in% names(mae_feature_data(mae, "rna")))

  annotation <- data.frame(
    chromosome = paste0("chr", seq_along(features)),
    biotype = rep(c("coding", "noncoding"), length.out = length(features)),
    row.names = rev(features),
    check.names = FALSE
  )
  annotation <- annotation[rev(rownames(annotation)), , drop = FALSE]
  with_annotation <- mae_add_feature_data(
    with_symbols,
    "rna",
    annotation
  )
  expect_identical(
    rownames(mae_feature_data(with_annotation, "rna")),
    features
  )
  expect_error(
    mae_add_feature_data(with_annotation, "rna", symbols, name = "symbol"),
    "already exist"
  )
  replaced <- mae_add_feature_data(
    with_annotation,
    "rna",
    stats::setNames(rep("replacement", length(features)), features),
    name = "symbol",
    overwrite = TRUE
  )
  expect_identical(
    mae_feature_data(replaced, "rna")$symbol,
    rep("replacement", length(features))
  )
  expect_error(
    mae_add_feature_data(mae, "rna", unname(symbols)),
    "name.*required"
  )
  expect_error(
    mae_add_feature_data(
      mae,
      "rna",
      symbols[-1L],
      name = "incomplete"
    ),
    "every feature"
  )
})

test_that("sample metadata bridges align assay samples without shadowing primary data", {
  mae <- .make_input_bridge_mae()
  sample_ids <- colnames(mae_pull_assay(mae, "rna", "counts"))
  risk <- stats::setNames(rev(seq_along(sample_ids) / 10), rev(sample_ids))

  with_risk <- mae_add_sample_data(
    mae,
    "rna",
    value = risk,
    name = "risk_score"
  )
  observed <- mae_samples(with_risk, "rna")
  expect_identical(rownames(observed), sample_ids)
  expect_equal(observed$risk_score, seq_along(sample_ids) / 10)
  expect_false("risk_score" %in% names(mae_samples(mae, "rna")))

  expect_error(
    mae_add_sample_data(
      mae,
      "rna",
      stats::setNames(rep("treated", length(sample_ids)), sample_ids),
      name = "condition",
      overwrite = TRUE
    ),
    "primary MAE metadata"
  )
  expect_error(
    mae_add_sample_data(with_risk, "rna", risk, name = "risk_score"),
    "already exist"
  )
  replaced <- mae_add_sample_data(
    with_risk,
    "rna",
    stats::setNames(rep(1, length(sample_ids)), sample_ids),
    name = "risk_score",
    overwrite = TRUE
  )
  expect_identical(mae_samples(replaced, "rna")$risk_score, rep(1, 4L))
})

test_that("sample metadata are keyed by raw assay columns under sampleMap", {
  primary_ids <- c("patient1", "patient2")
  assay_ids <- c("aliquotA", "aliquotB")
  counts <- matrix(
    seq_len(8L),
    nrow = 4L,
    dimnames = list(paste0("gene", 1:4), assay_ids)
  )
  leaf <- mae_create_experiment(
    counts,
    data.frame(row.names = assay_ids)
  )
  map <- data.frame(
    assay = "rna",
    primary = primary_ids,
    colname = assay_ids
  )
  mae <- mae_create(
    list(rna = leaf),
    data.frame(group = c("a", "b"), row.names = primary_ids),
    sample_map = map
  )
  purity <- stats::setNames(rev(c(0.7, 0.8)), rev(assay_ids))

  result <- mae_add_sample_data(mae, "rna", purity, name = "purity")
  raw_result <- MultiAssayExperiment::experiments(result)[["rna"]]

  expect_identical(rownames(SummarizedExperiment::colData(raw_result)), assay_ids)
  expect_equal(SummarizedExperiment::colData(raw_result)$purity, c(0.7, 0.8))
})

test_that("mae_add_experiment preserves maps, metadata, and the original object", {
  mae <- .make_input_bridge_mae()
  S4Vectors::metadata(mae)$study <- "simulated"
  MultiAssayExperiment::drops(mae) <- list(audit = "preserve")
  sample_ids <- rownames(MultiAssayExperiment::colData(mae))
  score <- matrix(
    seq_len(8L),
    nrow = 2L,
    dimnames = list(c("setA", "setB"), sample_ids)
  )

  expect_silent(
    with_scores <- mae_add_experiment(
      mae,
      score,
      name = "pathway",
      assay_name = "score"
    )
  )
  expect_identical(mae_experiments(with_scores), c("rna", "pathway"))
  expect_identical(mae_assays(with_scores, "pathway"), "score")
  expect_identical(mae_pull_assay(with_scores, "pathway", "score"), score)
  expect_identical(S4Vectors::metadata(with_scores)$study, "simulated")
  expect_identical(
    MultiAssayExperiment::drops(with_scores),
    list(audit = "preserve")
  )
  expect_identical(mae_experiments(mae), "rna")
  expect_error(
    mae_add_experiment(with_scores, score, name = "pathway"),
    "already exists"
  )
  replaced <- mae_add_experiment(
    with_scores,
    score * 2,
    name = "pathway",
    assay_name = "score",
    overwrite = TRUE
  )
  expect_identical(
    mae_pull_assay(replaced, "pathway", "score"),
    score * 2
  )

  invalid <- score
  colnames(invalid)[1L] <- "not_primary"
  expect_error(
    mae_add_experiment(mae, invalid, name = "invalid"),
    "must be primary sample IDs"
  )
})

test_that("mae_add_experiment can copy a source map or use an explicit map", {
  primary_ids <- c("patient1", "patient2")
  assay_ids <- c("aliquotA", "aliquotB")
  counts <- matrix(
    seq_len(8L),
    nrow = 4L,
    dimnames = list(paste0("gene", 1:4), assay_ids)
  )
  source <- mae_create_experiment(counts, data.frame(row.names = assay_ids))
  source_map <- data.frame(
    assay = "rna",
    primary = primary_ids,
    colname = assay_ids
  )
  mae <- mae_create(
    list(rna = source),
    data.frame(group = c("a", "b"), row.names = primary_ids),
    sample_map = source_map
  )
  protein <- matrix(
    seq_len(6L),
    nrow = 3L,
    dimnames = list(paste0("protein", 1:3), rev(assay_ids))
  )

  expect_silent(
    copied <- mae_add_experiment(
      mae,
      protein,
      name = "protein",
      source_experiment = "rna"
    )
  )
  copied_map <- as.data.frame(MultiAssayExperiment::sampleMap(copied))
  copied_map <- copied_map[as.character(copied_map$assay) == "protein", ]
  expect_setequal(copied_map$primary, primary_ids)
  expect_setequal(copied_map$colname, assay_ids)
  expect_identical(
    rownames(mae_samples(copied, "protein")),
    primary_ids
  )

  explicit <- mae_add_experiment(
    mae,
    protein,
    name = "protein",
    sample_map = data.frame(
      primary = rev(primary_ids),
      colname = rev(assay_ids)
    )
  )
  explicit_map <- as.data.frame(MultiAssayExperiment::sampleMap(explicit))
  explicit_map <- explicit_map[as.character(explicit_map$assay) == "protein", ]
  expect_identical(explicit_map$primary, rev(primary_ids))
  expect_identical(explicit_map$colname, rev(assay_ids))

  protein_se <- mae_create_experiment(
    protein,
    data.frame(platform = rep("panel", 2L), row.names = colnames(protein)),
    assay_name = "abundance"
  )
  from_se <- mae_add_experiment(
    mae,
    protein_se,
    name = "protein",
    source_experiment = "rna"
  )
  expect_identical(mae_assays(from_se, "protein"), "abundance")
})

test_that("adjust_covariates standardizes SVA, RUV, and matrix inputs", {
  samples <- paste0("sample", 1:4)
  sva_result <- list(
    sv = matrix(
      seq_len(8L) / 10,
      nrow = 4L,
      dimnames = list(samples, NULL)
    )
  )
  sva_covariates <- adjust_covariates(sva_result)
  expect_s3_class(sva_covariates, "data.frame")
  expect_identical(rownames(sva_covariates), samples)
  expect_identical(names(sva_covariates), c("SV1", "SV2"))
  expect_identical(
    adjust_covariates(sva_result, columns = "SV2"),
    sva_covariates[, "SV2", drop = FALSE]
  )

  ruv_result <- list(
    W = matrix(seq_len(8L) / 20, nrow = 4L),
    normalizedCounts = matrix(
      seq_len(24L),
      nrow = 6L,
      dimnames = list(paste0("gene", 1:6), samples)
    )
  )
  colnames(ruv_result$W) <- c("W_1", "W_2")
  ruv_covariates <- adjust_covariates(ruv_result, columns = 1L)
  expect_identical(rownames(ruv_covariates), samples)
  expect_identical(names(ruv_covariates), "W_1")

  matrix_input <- matrix(
    seq_len(8L),
    nrow = 4L,
    dimnames = list(samples, c("batch_score", "quality"))
  )
  expect_identical(
    adjust_covariates(matrix_input),
    as.data.frame(matrix_input)
  )
  expect_error(
    adjust_covariates(unname(matrix_input)),
    "real sample row names"
  )
  expect_error(
    adjust_covariates(list(sv = unname(sva_result$sv))),
    "real sample row names"
  )
})

test_that("adjust_covariates reads W columns from a SeqExpressionSet", {
  suppressWarnings(skip_if_not_installed("EDASeq"))
  suppressWarnings(skip_if_not_installed("Biobase"))
  samples <- paste0("sample", 1:4)
  counts <- matrix(
    seq_len(24L),
    nrow = 6L,
    dimnames = list(paste0("gene", 1:6), samples)
  )
  phenotype <- Biobase::AnnotatedDataFrame(data.frame(
    group = rep(c("a", "b"), each = 2L),
    W_1 = seq_len(4L) / 10,
    W_2 = seq_len(4L) / 20,
    row.names = samples,
    check.names = FALSE
  ))
  constructor <- getExportedValue("EDASeq", "newSeqExpressionSet")
  result <- constructor(counts, phenoData = phenotype)

  observed <- adjust_covariates(result)

  expect_identical(rownames(observed), samples)
  expect_identical(names(observed), c("W_1", "W_2"))
  expect_true(all(vapply(observed, is.numeric, logical(1))))
})

test_that("SVA covariates can be extracted and added back to an MAE", {
  suppressWarnings(skip_if_not_installed("sva"))
  set.seed(811L)
  sample_ids <- paste0("sample", seq_len(12L))
  feature_ids <- paste0("gene", seq_len(80L))
  condition <- factor(rep(c("control", "treated"), each = 6L))
  batch <- factor(rep(c("batch1", "batch2"), times = 6L))
  expression <- matrix(
    stats::rnorm(length(feature_ids) * length(sample_ids)),
    nrow = length(feature_ids),
    dimnames = list(feature_ids, sample_ids)
  )
  expression[seq_len(10L), condition == "treated"] <-
    expression[seq_len(10L), condition == "treated"] + 1
  samples <- data.frame(
    condition = condition,
    batch = batch,
    row.names = sample_ids
  )
  mae <- mae_from_matrix(
    expression,
    samples,
    experiment = "rna",
    assay = "log_expression"
  )

  result <- suppressMessages(adjust_sva(
    mae,
    "rna",
    assay = "log_expression",
    full = ~ condition + batch,
    null = ~ batch,
    n_surrogates = 1L
  ))
  covariates <- adjust_covariates(result)
  with_sva <- mae_add_sample_data(mae, "rna", covariates)

  expect_identical(rownames(covariates), sample_ids)
  expect_identical(names(covariates), "SV1")
  expect_identical(mae_samples(with_sva, "rna")$SV1, covariates$SV1)
})
