test_that("edgeR wrapper runs the complete quasi-likelihood workflow", {
  skip_if_not_installed("edgeR", minimum_version = "4.10.0")
  mae <- make_model_mae()

  fit <- de_edger(
    mae,
    "rna",
    formula = ~ condition,
    filter = FALSE
  )
  expect_s3_class(fit, "DGEGLM")
  expect_true(!is.null(fit$dispersion))

  tested <- de_edger(
    mae,
    "rna",
    formula = ~ condition,
    filter = FALSE,
    coef = "conditiontreated"
  )
  expect_s3_class(tested, "DGELRT")
})

test_that("DESeq2 wrapper builds and fits a name-aligned dataset", {
  skip_if_not_installed("DESeq2")
  mae <- make_model_mae()
  fit <- de_deseq2(
    mae,
    "rna",
    design = ~ condition,
    fitType = "mean",
    quiet = TRUE
  )

  expect_s4_class(fit, "DESeqDataSet")
  expect_identical(colnames(fit), paste0("sample", 1:8))
})

test_that("limma wrapper keeps expression features and sample alignment", {
  skip_if_not_installed("limma")
  mae <- make_model_mae()
  fit <- de_limma(mae, "rna", ~ condition, assay = "log_expression")

  expect_s3_class(fit, "MArrayLM")
  expect_identical(rownames(fit$coefficients), paste0("gene", 1:100))
})

test_that("GSVA uses the parameter-object API", {
  skip_if_not_installed("GSVA", minimum_version = "2.0.0")
  mae <- make_model_mae()
  sets <- list(
    first = paste0("gene", 1:20),
    second = paste0("gene", 21:40)
  )
  scores <- score_gsva(
    mae,
    "rna",
    gene_sets = sets,
    assay = "log_expression",
    kcdf = "Gaussian",
    min_size = 5L,
    verbose = FALSE
  )

  expect_identical(dim(scores), c(2L, 8L))
  expect_identical(colnames(scores), paste0("sample", 1:8))
})

test_that("decoupleR maps an explicit weight column", {
  skip_if_not_installed("decoupleR", minimum_version = "2.0.0")
  mae <- make_model_mae()
  network <- data.frame(
    source = rep(c("regulator1", "regulator2"), each = 10),
    target = paste0("gene", 1:20),
    weight = rep(c(-1, 1), 10),
    stringsAsFactors = FALSE
  )
  result <- activity_decouple(
    mae,
    "rna",
    network = network,
    assay = "log_expression",
    statistics = "ulm",
    method_args = list(NULL),
    consensus = FALSE,
    mor = "weight",
    min_size = 5L
  )

  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0L)
})

test_that("consensus clustering returns classes and diagnostics", {
  skip_if_not_installed("ConsensusClusterPlus")
  mae <- make_model_mae(n_features = 50L, n_samples = 8L)
  fit <- cluster_consensus(
    mae,
    "rna",
    assay = "log_expression",
    max_k = 2L,
    repetitions = 5L,
    item_fraction = 1,
    feature_fraction = 1,
    seed = 4L,
    title = tempdir(),
    plot = NULL
  )

  expect_true(is.list(fit))
  expect_length(cluster_consensus_classes(fit, 2L), 8L)
  diagnostics <- cluster_consensus_diagnostics(fit, title = tempdir(), plot = NULL)
  expect_true(is.list(diagnostics))
})

test_that("goseq accepts a named logical selection vector", {
  skip_if_not_installed("goseq")
  genes <- paste0("gene", seq_len(200))
  selected <- stats::setNames(seq_along(genes) <= 40L, genes)
  bias <- stats::setNames(seq_along(genes), genes)
  mapping <- data.frame(
    gene = genes,
    category = rep(paste0("set", 1:4), each = 50),
    stringsAsFactors = FALSE
  )

  result <- enrich_goseq(
    selected,
    genome = "custom",
    id = "custom",
    bias = bias,
    gene_to_category = mapping,
    categories = NULL,
    plot_fit = FALSE
  )
  expect_s3_class(result, "data.frame")
})

test_that("maSigPro rejects design rows that cannot align to samples", {
  skip_if_not_installed("maSigPro")
  mae <- make_model_mae(n_features = 20L, n_samples = 8L)
  design <- data.frame(
    Time = rep(c(0, 1), 4),
    Replicates = rep(1:4, each = 2),
    control = rep(c(1, 0), each = 4),
    treated = rep(c(0, 1), each = 4),
    row.names = paste0("wrong", 1:8)
  )

  expect_error(
    de_masigpro(
      mae,
      "rna",
      design,
      assay = "counts",
      counts = TRUE
    ),
    "sample"
  )
})

test_that("WGCNA quality filtering is automatic and auditable", {
  skip_if_not_installed("WGCNA")
  mae <- make_model_mae(n_features = 80L, n_samples = 12L)
  expression <- mae_pull_assay(mae, "rna", "log_expression")
  expression[1, ] <- 1
  expression[, 1] <- NA_real_
  samples <- mae_samples(mae, "rna")
  experiment <- mae_create_experiment(expression, samples, assay_name = "log")
  filtered_mae <- mae_create(list(rna = experiment), samples)

  expect_warning(
    result <- coexpr_wgcna(
      filtered_mae,
      "rna",
      "log",
      power = 2,
      min_module_size = 5,
      verbose = 0
    ),
    "removed"
  )
  expect_true("gene1" %in% result$bulkMAERemovedFeatures)
  expect_true("sample1" %in% result$bulkMAERemovedSamples)
  audit_names <- c(
    "bulkMAEQuality",
    "bulkMAERemovedSamples",
    "bulkMAERemovedFeatures",
    "bulkMAENetworkType"
  )
  expect_true(all(audit_names %in% names(result)))
})

test_that("module preservation uses a three-way feature intersection", {
  skip_if_not_installed("WGCNA")
  reference <- make_model_mae(n_features = 80L, n_samples = 12L, seed = 1L)
  test <- make_model_mae(n_features = 80L, n_samples = 12L, seed = 2L)
  reference_matrix <- mae_pull_assay(reference, "rna", "log_expression")
  test_matrix <- mae_pull_assay(test, "rna", "log_expression")
  reference_matrix[1, ] <- 1
  test_matrix[, 1] <- NA_real_
  reference_data <- mae_samples(reference, "rna")
  test_data <- mae_samples(test, "rna")
  reference <- mae_create(
    list(
      rna = mae_create_experiment(
        reference_matrix,
        reference_data,
        assay_name = "log"
      )
    ),
    reference_data
  )
  test <- mae_create(
    list(
      rna = mae_create_experiment(
        test_matrix,
        test_data,
        assay_name = "log"
      )
    ),
    test_data
  )
  module_colors <- stats::setNames(
    rep(c("blue", "brown"), each = 30L),
    paste0("gene", 1:60)
  )

  expect_warning(
    result <- coexpr_preservation(
      reference,
      "rna",
      test,
      "rna",
      module_colors,
      reference_assay = "log",
      test_assay = "log",
      permutations = 2L,
      verbose = 0
    ),
    "removed"
  )
  audit <- result$bulkMAEInputQuality
  expect_true("sample1" %in% audit$removedTestSamples)
  expect_false("gene1" %in% audit$analyzedFeatures)
  expect_true(all(paste0("gene", 61:80) %in% audit$excludedReferenceFeatures))
  expect_identical(audit$networkType, "signed")
})
